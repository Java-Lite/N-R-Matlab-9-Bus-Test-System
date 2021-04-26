clear all
clc
%% Bus Datas
% % %     Bus Bus  Voltage  Angle  ---Generator------Load---
% % %     No  code   Mag.   Degree    MW    Mvar   MW    Mvar 
busdata=  [1   1    1.040     0        0     0      0      0;
           2   2    1.025     0       163    0      0      0;
           3   2    1.025     0        85    0      0      0;
           4   3    1.000     0        0     0      0      0;
           5   3    1.000     0        0     0     125    50;
           6   3    1.000     0        0     0      90    30;
           7   3    1.000     0        0     0      0      0;
           8   3    1.000     0        0     0     100    35;
           9   3    1.000     0        0     0      0      0];
%% Line Datas                                  
% %      Bus  bus     R        X       1/2 B    
% %       fb  tb     p.u.     p.u.     p.u.       T
linedata=[1   4    0.00000   0.05760  0.00000     1;
          2   7    0.00000   0.06250  0.00000     1;
          3   9    0.00000   0.05860  0.00000     1;
          4   5    0.01000   0.08500  0.17600/2   1;
          4   6    0.01700   0.09200  0.15800/2   1;
          5   7    0.03200   0.16100  0.30600/2   1;
          6   9    0.03900   0.17000  0.35800/2   1;
          7   8    0.00850   0.07200  0.14900/2   1;
          8   9    0.01190   0.10080  0.20900/2   1];
%% Data arranged for Linedata in the different vector 
fb=linedata(:,1);tb=linedata(:,2);
r=linedata(:,3);x=linedata(:,4);
b=linedata(:,5);a=linedata(:,6);
z=r+1i*x; 						% Impedance of Branch 
y=1./z;b=1i*b;  				% Admittance of Branch 
nl=length(fb);					% No of Branch 
No_of_Bus=max(max(fb),max(tb)); % No of Bus 
%% Formation of YBus matrix 
Y=zeros(No_of_Bus,No_of_Bus);				
for k=1:nl
    Y(fb(k),tb(k))=Y(fb(k),tb(k))-y(k)/a(k);
    Y(tb(k),fb(k))=Y(fb(k),tb(k));
end
for m=1:No_of_Bus                           
    for n=1:nl
        if fb(n)==m
            Y(m,m)=Y(m,m)+y(n)/a(n)^2+b(n);
        elseif tb(n)==m
            Y(m,m)=Y(m,m)+y(n)+b(n);
        end
    end
end
AP=real(Y);
RP=imag(Y);
G=abs(Y);
B=angle(Y)/pi*180; 
%% Initializations Bus Datas
BMva=100;
busNo=busdata(:,1);type=busdata(:,2);V=busdata(:,3);del=busdata(:,4);
Pg=busdata(:,5)/BMva;Qg=busdata(:,6)/BMva;Pl=busdata(:,7)/BMva;
Ql=busdata(:,8)/BMva;PV_Bus=find(type==2|type==1);
PQ_Bus=find(type==3);% type1(Slack),type2(PV_Bus Bus),type3(PQ_Bus Bus )
No_of_PQ_Bus=length(PQ_Bus);No_of_PV_Bus=length(PV_Bus);
Active_Power=Pg-Pl;Reactive_Power=Qg-Ql; 
Iter=1;Tol=1; % Iterantion And tolerance 
%% Newton Raphson Load Flow 
while Tol>4e-6
P = zeros(nl,1);
Q = zeros(nl,1);              % Calculate P and Q
      for i = 1:nl
          for k = 1:nl
              P(i)=P(i)+V(i)*V(k)*G(i,k)*cosd(B(i,k)-del(i)+del(k));
              Q(i)=Q(i)-V(i)*V(k)*G(i,k)*sind(B(i,k)-del(i)+del(k));
          end
      end
Psch=Active_Power-P;
Qsch=Reactive_Power-Q;
dP=Psch(2:No_of_Bus);         % Looking for Power Residuals of (P) in p.u
k=1;
dQ=zeros(No_of_PQ_Bus,1);
    for i=1:No_of_Bus
        if type(i)==3
            dQ(k,1)=Qsch(i);  % Looking for Power Residuals of (Q) in p.u
            k=k+1;
        end
    end
 M=[dP;dQ];                   % Delta Matrix in p.u
%% Formation Of J1
 J1 = zeros(No_of_Bus-1,No_of_Bus-1);
 for i = 1:(No_of_Bus-1)
     m = i+1;
     for k = 1:(No_of_Bus-1)
         n = k+1;
         if n == m
            for n=1:No_of_Bus
                J1(i,k) = J1(i,k) + V(m)*V(n)*G(m,n)*sind(B(m,n)-del(m)+del(n));
            end
                J1(i,k) = J1(i,k)-V(m)^2*RP(m,m);
         else
            J1(i,k) = J1(i,k)-V(n)*V(m)*G(n,m)*sind(B(n,m)-del(n)+del(m));
         end
      end
 end
%% Formation Of J2
 J2=zeros(No_of_Bus-1,No_of_PQ_Bus);
 for i=1:No_of_Bus-1
     m=i+1;
     for j=1:No_of_PQ_Bus
         n=PQ_Bus(j);
         if m==n
             for n=1:No_of_Bus
                 J2(i,j) = J2(i,j) + V(m)*V(n)*G(m,n)*cosd(B(m,n)-del(m)+del(n));
             end
                 J2(i,j)=J2(i,j)+V(m)*AP(m,m);
         else
                 J2(i,j) = V(m)*V(n)*G(m,n)*cosd(B(m,n)-del(m)+del(n));
         end
     end
 end
%% Formation Of J3
 J3 = zeros(No_of_PQ_Bus,No_of_Bus-1);
 for i = 1:No_of_PQ_Bus
     m = PQ_Bus(i);
     for j = 1:No_of_Bus-1
         n = j+1;
         if n == m
              for n=1:No_of_Bus
                 J3(i,j) = J3(i,j) + V(m)*V(n)*G(m,n)*cosd(B(m,n)-del(m)+del(n));
              end
                 J3(i,j)=J3(i,j)-V(m)^2*AP(m,m);
         else
                 J3(i,j) = -V(n)*V(m)*G(n,m)*cosd(B(n,m)-del(n)+del(m));
         end
     end
 end
%% Formation Of J4
J4=zeros(No_of_PQ_Bus,No_of_PQ_Bus);
for i=1:No_of_PQ_Bus
    m=PQ_Bus(i);
    for j=1:No_of_PQ_Bus
        n=PQ_Bus(j);
        if m==n
            for n=1:No_of_Bus
                J4(i,j) = J4(i,j) - V(n)*G(m,n)*sind(B(m,n)-del(m)+del(n));
            end
                J4(i,j)=J4(i,j)-V(m)*RP(m,m);
        else
                J4(i,j) = - V(n)*G(m,n)*sind(B(m,n)-del(m)+del(n));
        end
    end
end
J=[J1 J2;J3 J4];
X=inv(J)*M;
dTh=X(1:No_of_Bus-1);                   % Change in angle 
dV=X(No_of_Bus:end);                    % change in volatge mag 
del(2:No_of_Bus)=del(2:No_of_Bus)+dTh;  % Voltage angle update 

% Voltage mag update 
k=1;
for n=2:No_of_Bus
   if type(n)==3
      V(n)=V(n)+dV(k);
      k=k+1;
   end
end
    Convergence_curve=V;
    Iter=Iter+1;
    Tol=max(abs(M));
end
% Conversion del from radian to degree
deli=del.*3.14/180;

%% Arus saluran from bus to bus
VI=V+1i.*deli;
for n=1:nl
    If(n)=((VI(fb(n))-VI(tb(n))).*y(n));
    Ir(n)=((VI(tb(n))-VI(fb(n))).*y(n));
end
   If=real(If)+1i*(imag(If).*-1);
   Ir=real(Ir)+1i*(imag(Ir).*-1);

%% Aliran daya dan Rugi Daya
for n=1:nl
    Sf(n)=(VI(fb(n))).*(If(n));
    Sr(n)=(VI(tb(n))).*(Ir(n));
end
   V;
   del;  
   P=P*BMva;
   Q=Q*BMva;
   Qr=sum(abs(Q));
   S=sqrt((P.^2) + (Q.^2));
   Pf=P./S;
   Pf_mean=mean(abs(Pf));
   If;
   Ir;
   Sf=Sf.*BMva; 
   Sr=Sr.*BMva; 
   SSS=Sf+Sr;
   SS=sum(SSS);
   V=V.*BMva;  
%% Looking for Ploss
for n=1:nl % No of Branch
    Pl(n)=((P(tb(n)).^2+Q(tb(n)).^2)./V(tb(n)).^2).*r(n);
    Ql(n)=((P(tb(n)).^2+Q(tb(n)).^2)./V(tb(n)).^2).*x(n);
end
for n=1:nl % No of Branch
    Pl(n)=(Pl(n))./Q(tb(n));
    Ql(n)=(Ql(n))./P(tb(n));
end
V=V./BMva;
% Looking for Ploss
for n=1:nl % No of Branch
    Pll(n)=(2.*Q(tb(n)).*r(n))./V(tb(n)).^2;
    Qll(n)=(2.*Q(tb(n)).*x(n))./V(tb(n)).^2;
end
Pll=Pll';
Qll=Qll';
Vnorm=V./0.95;
 %% Load Flow Solution 
disp('----------------------------------------');
disp('  Newton Raphson Loadflow Solution    ');
disp('----------------------------------------');
disp(' |Bus |   |Voltage|    |Angle |   |Power |  |Active|  |Reactive|');
disp(' | No |   |pu     |    |Degree|   |Factor|  |Power |  |Power   |');
disp('-------------------------------------------------------------------');
for m=1:No_of_Bus 
    fprintf(' %4g   ' ,m);
    fprintf(' %8.3f    ' ,V(m));
    fprintf(' %8.3f  ' ,del(m));
    fprintf(' %8.3f  ' ,Pf(m));
    fprintf(' %8.3f  ' ,P(m));
    fprintf(' %8.3f  ' ,Q(m));
    fprintf('\n');
end
disp('-------------------------------------------------------------------');
display(['System Losses ', num2str(SS)]);
fprintf( 'Number Of Ieration %4g \n',Iter); 
 %% Loss Sensitivity Factors 
disp('----------------------------------------');
disp('  Loss Sensitivity Factors    ');
disp('----------------------------------------');
disp(' |Voltages | |Active Power| |Reactive Power|');
disp(' |Norm     | |Losses      | |Losses        |');
disp(' |(Bus)    | |(Line)      | |(Line)        |');
disp('----------------------------------------');
for m=1:No_of_Bus 
    fprintf(' %9.4f   ' ,Vnorm(m));
    fprintf(' %9.4f    ' ,Pll(m));
    fprintf(' %9.4f  ' ,Qll(m));
    fprintf('\n');
end
disp('----------------------------------------');
disp('### Created by Irfan Rizali ###');
%% Draw objective space
subplot(1,1,1);
semilogy(Convergence_curve,'Color','r')
title('All Voltages')
xlabel('Bus');
ylabel('Voltages of the Bus');
axis tight
grid on
box on
legend('Volt')