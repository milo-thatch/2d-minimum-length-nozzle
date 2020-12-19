%Ludovico Fossà 12/2020
%Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
function nu=prandtl_meyer(mach,gamma)
    nu=sqrt((gamma+1)./(gamma-1)).*atan(sqrt((gamma-1)./(gamma+1).*(mach.^2-1)))-atan(sqrt(mach.^2-1));
end