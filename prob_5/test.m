clear all
close all

load('points.dat');


scatter(points(:,1),points(:,2),50,points(:,3),'fill')
colorbar