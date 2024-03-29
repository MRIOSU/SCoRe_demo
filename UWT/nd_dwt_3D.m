%Calculates the nondecimated Wavelet transform using fast convlolution and
%   periodic boundry conditions.  The Filters are pre-computed and stored
%   in f_dec and f_rec to reduce computation time when performing Multiple
%   forward/backward transforms
%
%Methods:
%   ud_dwt3D:   Constructor
%               Inputs: wname - Wavelet Filters to Use i.e. db1,db2,etc.
%                        Either a string 'db1' or a cell {'db4','db1'} where
%                        the first element in the string is filter for the 
%                        spatial domain and the second is the filter for the 
%                        time domain
%
%                       sizes - size of the 3D object [n1,n2,n3]
%
%                       preserve_l2_norm - An optional third input.  If set
%                        TRUE, the l2 norm in the wavelet domain will be
%                        equal to the l2 norm in the signal domain
%
%   dec:        Multilevel Decomposition
%               Inputs: x - Image domain signal for decomposition
%
%                       levels - Number of decomposition Levels
%
%               Outputs: y - Multilevel non-decimated DWT coefficients in a
%                        4D array where the data is arranged [n1,n2,n3,bands]
%                        The bands are orginized as follows.  Let "HVT"
%                        represent the Horizontal, Vertical, and Temporal
%                        Bands.  The coefficients are ordered as follows,
%                        "LLL", "HLL", "LHL", "HHL", "LLH", "HLH", "LHH",
%                        "HHH" where "H" denotes the high frequency filter
%                        and L represents the low frequency filter.
%                        Successive levels of decomposition are stacked
%                        such that the highest "LLL" is in [n1,n2,n3,1]
%
%   rec:        Multilevel Reconstruction
%               Inputs: x - Wavelet coefficients in a 4D array size 
%                       [n1,n2,n3,bands].
%
%               Outputs: y - Reconstructed 3D array.%   
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Last update:  2/5/2015
%**************************************************************************

classdef nd_dwt_3D
    properties
        f_dec;          % Decomposition Filters
        sizes;          % Size of the 3D Image
        f_size;         % Length of the filters
        wname;          % Wavelet used
        pres_l2_norm;   % Binary indicator to preserver l2 norm of coefficients
    end
    
    %% Public Methods
    methods
        % Constructor Object
        function obj = nd_dwt_3D(wname,sizes,varargin)
            % Set Image size
            if length(sizes) ~=3
                error('The sizes vector must be length 3');
            else
                obj.sizes = sizes;
            end
            
            if ischar(wname)
                obj.wname = {wname,wname,wname};
            elseif iscell(wname)
                if length(wname) ==3
                    obj.wname = wname;
                else
                    error(['You must specify three filter names in a cell array'...
                            ,'of length 3, or a single string for the same'...
                            ,' filter to be used in all dimensions']);
                end
            end
            
            if isempty(varargin)
                obj.pres_l2_norm = 0;
            else
                obj.pres_l2_norm = varargin{1};
            end

            % Get the Filter Coefficients
            [obj.f_dec,obj.f_size] = obj.get_filters(obj.wname);
            
        end
        
        % Multilevel Undecimated Wavelet Decomposition
        function y = dec(obj,x,level)
            % Fourier Transform of Signal
            x = fftn(x);
            
            % Preallocate
            y = zeros([obj.sizes, 8+7*(level-1)]);
            
            % Calculate Mutlilevel Wavelet decomposition
            for ind = 1:level
                % First Level
                if ind ==1
                    y = level_1_dec(obj,x);
                % Succssive Levels
                else
                    y = cat(4,level_1_dec(obj,fftn(squeeze(y(:,:,:,1)))), y(:,:,:,2:end));
                end
            end
                      
        end
        
        % Multilevel Undecimated Wavelet Reconstruction
        function y = rec(obj,x,level) % Rizwan: Insert 'level'
            
            bands = 7*level+1; % Rizwan: Add these two lines
            x = reshape(x,[obj.sizes, bands]);
                       
            % Fourier Transform of Signal
            x = fft(fft(fft(x,[],1),[],2),[],3);
            
            % Reconstruct from Multiple Levels
            for ind = 1:level
                % First Level
                if ind ==1
                    y = level_1_rec(obj,x);
                    if ~obj.pres_l2_norm
                        y = y/8;
                    end
                % Succssive Levels
                else
                    y = fftn(y);
                    y = level_1_rec(obj,cat(4,y,x(:,:,:,9+(ind-2)*7:15+(ind-2)*7)));
                    if ~obj.pres_l2_norm
                        y = y/8;
                    end
                end
            end 
            
        end
    end
    
    %% Private Methods
    methods (Access = protected,Hidden = true)
        
        % Returns the Filters and 
        function [f_dec,f_size] = get_filters(obj,wname)
        % Decomposition Filters
        
            % Get Filters for Spatial Domain
            [LO_D,HI_D] = wave_filters(wname{1});
            [LO_D2,HI_D2] = wave_filters(wname{2});
            [LO_D3,HI_D3] = wave_filters(wname{3});
            
            % Find the filter size
            f_size.s1 = length(LO_D);
            f_size.s2 = length(LO_D2);
            f_size.s3 = length(LO_D3);
            
            % Dimension Check
            if f_size.s1 > obj.sizes(1)
                error(['First Dimension of Data is shorter than the wavelet'...
                    ,' filter being used']);
            elseif f_size.s2 > obj.sizes(2)
                error(['Second Dimension of Data is shorter than the wavelet'...
                    ,' filter being used']);
            elseif f_size.s3 > obj.sizes(3)
                error(['Third Dimension of Data is shorter than the wavelet'...
                    ,' filter being used']);
            end
            
            % Get the 2D Filters by taking outer products
            dec_LL = LO_D.'*LO_D2;
            dec_HL = HI_D.'*LO_D2;
            dec_LH = LO_D.'*HI_D2;
            dec_HH = HI_D.'*HI_D2;

            % Take the Outerproducts for the third dimension
            for ind = 1:size(dec_LL,2)
                f_dec.LLL(:,ind,:) = dec_LL(:,ind)*LO_D3;
                f_dec.HLL(:,ind,:) = dec_HL(:,ind)*LO_D3;
                f_dec.LHL(:,ind,:) = dec_LH(:,ind)*LO_D3;
                f_dec.HHL(:,ind,:) = dec_HH(:,ind)*LO_D3;
                f_dec.LLH(:,ind,:) = dec_LL(:,ind)*HI_D3;
                f_dec.HLH(:,ind,:) = dec_HL(:,ind)*HI_D3;
                f_dec.LHH(:,ind,:) = dec_LH(:,ind)*HI_D3;
                f_dec.HHH(:,ind,:) = dec_HH(:,ind)*HI_D3;
            end

            % Add a circularshift of half the filter length to the 
            % reconstruction filters by adding phase to them
            phase1 = exp(1j*2*pi*f_size.s1/2*linspace(0,1-1/obj.sizes(1),obj.sizes(1)));
            phase2 = exp(1j*2*pi*f_size.s2/2*linspace(0,1-1/obj.sizes(2),obj.sizes(2)));
            phase3 = exp(1j*2*pi*f_size.s3/2*linspace(0,1-1/obj.sizes(3),obj.sizes(3)));

            % 2D Phase
            shift_t = phase1.'*phase2;
            
            % 3D Phase
            shift = zeros([obj.sizes(1),obj.sizes(2),obj.sizes(3)]);
            for ind = 1:size(shift_t,2)
               shift(:,ind,:) = shift_t(:,ind)*phase3; 
            end
            
            % Take the Fourier Transform of the Kernels for Fast
            % Convolution
            if obj.pres_l2_norm 
                scale = 1/sqrt(8);
            else
                scale = 1;
            end
            f_dec.LLL = scale*shift.*fftn(f_dec.LLL,[obj.sizes(1),obj.sizes(2),obj.sizes(3)]);
            f_dec.HLL = scale*shift.*fftn(f_dec.HLL,[obj.sizes(1),obj.sizes(2),obj.sizes(3)]);
            f_dec.LHL = scale*shift.*fftn(f_dec.LHL,[obj.sizes(1),obj.sizes(2),obj.sizes(3)]);
            f_dec.HHL = scale*shift.*fftn(f_dec.HHL,[obj.sizes(1),obj.sizes(2),obj.sizes(3)]);
            f_dec.LLH = scale*shift.*fftn(f_dec.LLH,[obj.sizes(1),obj.sizes(2),obj.sizes(3)]);
            f_dec.HLH = scale*shift.*fftn(f_dec.HLH,[obj.sizes(1),obj.sizes(2),obj.sizes(3)]);
            f_dec.LHH = scale*shift.*fftn(f_dec.LHH,[obj.sizes(1),obj.sizes(2),obj.sizes(3)]);
            f_dec.HHH = scale*shift.*fftn(f_dec.HHH,[obj.sizes(1),obj.sizes(2),obj.sizes(3)]);
        end
        
        % Single Level Redundant Wavelet Decomposition
        function y = level_1_dec(obj,x_f)
            % Preallocate
            y = zeros([obj.sizes,8]);
            
            % Calculate Wavelet Coefficents Using Fast Convolution
            y(:,:,:,1) = ifftn(x_f.*obj.f_dec.LLL); 
            y(:,:,:,2) = ifftn(x_f.*obj.f_dec.HLL);
            y(:,:,:,3) = ifftn(x_f.*obj.f_dec.LHL);
            y(:,:,:,4) = ifftn(x_f.*obj.f_dec.HHL);
            y(:,:,:,5) = ifftn(x_f.*obj.f_dec.LLH);
            y(:,:,:,6) = ifftn(x_f.*obj.f_dec.HLH);
            y(:,:,:,7) = ifftn(x_f.*obj.f_dec.LHH);
            y(:,:,:,8) = ifftn(x_f.*obj.f_dec.HHH);       
        end
        
        % Single Level Redundant Wavelet Reconstruction
        function y = level_1_rec(obj,x_f)
            
            % Reconstruct the 3D array using Fast Convolution
            y = (squeeze(x_f(:,:,:,1)).*conj(obj.f_dec.LLL));
            y = y + (squeeze(x_f(:,:,:,2)).*conj(obj.f_dec.HLL));
            y = y + (squeeze(x_f(:,:,:,3)).*conj(obj.f_dec.LHL));
            y = y + (squeeze(x_f(:,:,:,4)).*conj(obj.f_dec.HHL));
            y = y + (squeeze(x_f(:,:,:,5)).*conj(obj.f_dec.LLH));
            y = y + (squeeze(x_f(:,:,:,6)).*conj(obj.f_dec.HLH));
            y = y + (squeeze(x_f(:,:,:,7)).*conj(obj.f_dec.LHH));
            y = y + (squeeze(x_f(:,:,:,8)).*conj(obj.f_dec.HHH));
            
            y = ifftn(y);
            
        end
    end
end

