%% Triangle Threshold Function
function [thresh]=triangle_th_david(hist, num_bins)
%        Method: 1) A straight line is drawn from the maximum to the end of the histogram, to form a triangular-like shape 
%                2) The threshold is selected at the point of the histogram that maximizes the perpendicular distance from 
%                   the histogram to the straight line
%
%        Inputs:    histogram -   histogram of the gray level image
%                    num_bins -   number of bins (e.g. gray levels is 256)
%       Outputs:       thresh -   threshold value in the range [0 1];
%         
%    Contraints: This method gives suitable results in many cases, but tends to be sensitive to parameters such as the 
%                statistical fluctuations of the histogram, and the position of the endpoint (highest gray level).
%
%    References: Zack (Zack GW, Rogers WE, Latt SA (1977), "Automatic measurement of sister chromatid exchange frequency",
%                J. Histochem. Cytochem. 25 (7): 741,53, )


thresh = 0;

%   Find maximum of histogram and its location along the x axis
[maxval, maxloc] = max(hist);
   
%   Find location of the left end of the histogram.
nonzerolocs = find(hist>0);
minloc = nonzerolocs(1);
  
% Calculate slope of hypotenuse of triangle - rise/run  
slope = maxval/(maxloc-minloc); 
    
% Use triangle geometry to calculate all possible distances d perpendicular from the hypotenuse to the histogram
hist = hist';
x1 = 0:(maxloc - minloc); % start x1 from 0 with the length of maxloc - minloc
y1 = hist(x1 + minloc);
beta = y1 + x1/slope;
x2 = beta / (slope + (1/slope));
y2 = slope * x2;

d = ((y2 - y1).^2 + (x2-x1).^2).^0.5; % Pythagoras Theorem: c = sqrt(a^2 + b^2)

%  Threshold is the location that maximizes distance d     
thresh = find(max(d)==d);
    
%   output of the program
thresh = minloc + mean(thresh);
thresh = thresh / num_bins;
    
end    
   
