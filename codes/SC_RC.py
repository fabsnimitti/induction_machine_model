from sympy import *
import math
import matplotlib.pyplot as plt

pi=math.pi
N=2 
Rs=.8467
Rr=.5176
Lm=66.0391/(2*pi*60)
Lls=2.1750/(2*pi*60)
Llr=2.5067/(2*pi*60)
Ls = Lls+Lm
Lr = Llr+Lm
J=0.0698
Bn = 0.015
P=2

sigma=1-(Lm**2)/(Ls*Lr)
neta=Rr/Lr
beta=Lm/(sigma*Ls*Lr)
gama=Rs/(sigma*Ls)+beta*neta*Lm
Lm_=Lm**2/Lr
eta=neta

u=[]

Ts=1/10000
tsim=3

counter=0
while(counter<=tsim):
    u.append(Matrix([[380*math.sin(2*pi*60*counter)],[380*math.cos(2*pi*60*counter)]]))

    counter=counter+Ts


B=Matrix([[1/(sigma*Ls) ,0],[0 ,1/(sigma*Ls)],[-beta, 0],[0, -beta]]);

x=[]
x.append(Matrix([[0],[0],[0],[0]]))
wr=[]
wr.append(0)
time=0
wa=[]
counter=0
pnt=[]
while(time<tsim):
    pnt.append(time)
    wa.append([0])
    A = Matrix([[Lr*Rs/(Lm**2-Ls*Lr), wr[counter]*Lm**2/(Lm**2-Lr*Ls), -Lm*Rr/(Lm**2-Lr*Ls), Lm*Lr*wr[counter]/(Lm**2-Lr*Ls)],
            [-wr[counter]*Lm**2/(Lm**2-Lr*Ls) ,Lr*Rs/(Lm**2-Ls*Lr) ,-wr[counter]*Lm*Lr/(Lm**2-Ls*Lr), -Lm*Rr/(Lm**2-Lr*Ls)],
            [ -Lm*Rs/(Lm**2-Lr*Ls), -wr[counter]*Lm*Ls/(Lm**2-Ls*Lr),Ls*Rr/(Lm**2-Lr*Ls),-wr[counter]*Lr*Ls/(Lm**2-Lr*Ls)],
            [wr[counter]*Lm*Ls/(Lm**2-Lr*Ls),-Lm*Rs/(Lm**2-Lr*Ls), wr[counter]*Lr*Ls/(Lm**2-Lr*Ls) ,Ls*Rr/(Lm**2-Lr*Ls)]])
    
    x.append((eye(4)+Ts*A)*x[counter] + Ts*B*u[counter])
    Te=(3/2)*(P/2)*(Lm/Lr)*(x[counter][0]*x[counter][3]-x[counter][1]*x[counter][2] );
    wr.append( (1-Ts*Bn/J)*wr[counter] + Ts/J*Te);
    time=time+Ts
    print(time)
    counter=counter+1

iqs=[]
ids=[]
iqr=[]
idr=[]
for j in x:
    iqs.append(j[0])
    ids.append(j[1])
    iqr.append(j[2])
    idr.append(j[3])
    
time_to_plot=[10000,12000]
plt.plot(iqs[time_to_plot[0]:time_to_plot[1]], 'g')
plt.legend(["Iqs", "Velocidade"], loc=0)
plt.show()



plt.plot(ids[time_to_plot[0]:time_to_plot[1]], 'b')
plt.legend(["Ids", "Velocidade"], loc=0)
plt.show()




plt.plot(iqr[time_to_plot[0]:time_to_plot[1]], 'black')
plt.legend(["Iqr", "Velocidade"], loc=0)
plt.show()



plt.plot(idr[time_to_plot[0]:time_to_plot[1]], 'violet')
plt.legend(["Idr", "Velocidade"], loc=0)
plt.show()



