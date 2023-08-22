function [segvector,segCM] = bodyseg(proxX,proxZ,distalX,distalZ,CMref)

% Body segment function
% by Rodrigo Bini - bini.rodrigo@gmail.com
% date Oct 2009
% -------------------------------------------------------------------------
% Calculate body segment and CM based on two spatial variables (distal and
% proximal) and CM reference. Exports the vector with the body segment size   
% and spatial coordinates of the centre of mass.
% =========================================================================

CMref = CMref/100;

segvector = sqrt(((proxZ-distalZ).^2)+((proxX-distalX).^2));
segCM = [proxX-((proxX-distalX).*CMref),proxZ-(((proxZ-distalZ).*CMref))];
