clear all

load('points.dat');

n = length(points);
step=floor(max(1,n/1000));

% scatter(points(1:step:end,1),points(1:step:end,2),20,points(1:step:end,3)/max(points(:,3)),'fill')
scatter(points(1:step:end,1),points(1:step:end,2),20,points(1:step:end,5),'fill')
colorbar


load('boxes.dat');

n_leaves = length(boxes); 

for i=1:n_leaves
    rectangle('Position',[boxes(i,2),boxes(i,3),boxes(i,4),boxes(i,5)]);
end


