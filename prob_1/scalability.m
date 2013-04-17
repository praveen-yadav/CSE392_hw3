clear all

time = [0.24620, 0.196327,  0.124098, 0.0663806,  0.0281576];
n_th = [2 4 8 16 32];

loglog(n_th, time, '-o')
xlabel('number of threads');
ylabel('time [s]');

title('run on stampede');
saveas(gca,'scalability_data.jpg');