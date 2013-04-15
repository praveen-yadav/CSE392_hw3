clear all

points(:,3) = points(:,3)/max(points(:,3));

scatter(points(:,1),points(:,2),50,points(:,3),'fill')
colorbar