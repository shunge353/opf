function mpc = case3
%CASE3

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	30.99	0	0	1	1	0	230	1	1.1	0.9;
	2	1	0	105.35	0	0	1	1	0	230	1	1.1	0.9;
	3	1	150	123.94	0	0	1	1	0	230	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	318	0	100	-100	1.02	100	1	200	0	0	0	0	0	0	0	0	0	0	0	0;
	2	0	0	100	-100	1	100	1	200	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.01008	0.1	0.1025	1000	250	250	0	0	1	-360	360;
	1	3	0.00744	0.3	0.0775	1000	250	250	0	0	1	-360	360;
	2	3	0.00744	0.1	0.0775	1000	250	250	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	2	60	0;
	2	0	0	2	120	0;
];
