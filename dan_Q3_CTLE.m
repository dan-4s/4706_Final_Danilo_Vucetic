i = 2*pi*10^10;
gm = 0.04; 
rd = 176; 
rs = 200; 
cs = 0.795775*10^-12; 
cp = 1/(i*rd);
a = gm/(cp);
b = (1+(gm*rs/2))/(rs*cs);
c = 1/(rd*cp);
d = 1/(rs*cs);

sys = tf([a, a*d],[1, b+c, b*c]);
bode(sys);
grid on;


