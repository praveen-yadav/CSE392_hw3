clear all
close all

load('points.dat');

n = length(points);
step=max(1,n/1000);

scatter(points(1:step:end,1),points(1:step:end,2),50,points(1:step:end,3)/max(points(:,3)),'fill')
colorbar