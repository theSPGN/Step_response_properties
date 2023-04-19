clc
clear all
close all
%Step could not start in (0,0) following lines provides correct definition
%of c0 and c1 that are next 2 max differences from final value
L= [1 5 5];
M= [1 1.65 5 6.5 2];
Percent_of_c0=0.02;
precision=100;
sys=tf(L,M);
[Pole,Damping,Frequency]=damp(sys);
dzeta = min(Damping);
omega = max(Frequency(Damping==dzeta));
tend = ceil( -log(Percent_of_c0) / abs(dzeta*omega) ); 
t = 1/precision:1/precision:tend;
[y,x] = lsim(L,M,ones(tend*precision,1),t);
[miny,maxy] = bounds(y-y(length(y)));
if abs(miny)<abs(maxy)
% if max> abs(min) then c0 max , c1 min for y that begins from
% value y(tmax:,1)
NewBegginingValue=maxy+y(length(y));
ynew = y(find(y==NewBegginingValue):length(y),1);
[miny,~] = bounds(ynew-y(length(y)));
end

if abs(miny)>abs(maxy)
% if abs(min)> abs(max) then c0 min, c1 max for y zaczthat begins from
% value y(tmin:,1)
NewBegginingValue=miny+y(length(y));
ynew = y(find(y==NewBegginingValue):length(y),1);
[~,maxy] = bounds(ynew-y(length(y)));
end




c0 = max(abs(miny),abs(maxy)); 
c1 = min(abs(miny),abs(maxy));
SettlingTime = 100*c1/c0;
delta_r= Percent_of_c0*abs(c0); % from 0.02 to 0.05 depends of assumptions 
wartosc_zadana = 1;
es=abs(y(length(y))-wartosc_zadana); % steady-state error

tr=0;
for i = t
index = int64(i*precision);
actual_value = y(index);
end_value = y(length(y));
r1 = actual_value - (end_value + delta_r); %finding last time when step
r2 = actual_value - (end_value - delta_r); %response cut the line
if r1 >= 0 ||  r2 <=0
tr = i;
end
end
plot(t,y,'r-');
hold on;
plot(t,y(length(y))*ones(tend*precision,1),'g-');

dr1=(y(length(y))+delta_r)*ones(tend*precision,1);
dr2=(y(length(y))-delta_r)*ones(tend*precision,1);
plot(t,[dr1,dr2],'k--');
xlim([0, tend])
ylim([0, max(ynew)*1.1])
xlabel('Time [s]')
ylabel('Step response')
legend('Response','Steady-state error','Precision*c0')
grid minor;

%properties of our object:
disp('es')
disp(es)
disp('c0')
disp(c0)
disp('c1')
disp(c1)
disp('Overshoot')
disp("    "+string(SettlingTime)+ " %")

%delta_r1 delta_r2 is defined % of c0
%tr - settling time
disp('delta_r')
disp(delta_r)
disp('Settling time')
disp(tr)

plot([0,tend],[max(ynew),max(ynew)],'Color',[0,0,0,0.25],'HandleVisibility','off')
plot([tr,tr],[0,max(ynew)*1.1],'Color',[0,0,0,0.25],'HandleVisibility','off')
text(tr,0.1*max(ynew),'Settling time','Color',[0.8,0.8,0.8])