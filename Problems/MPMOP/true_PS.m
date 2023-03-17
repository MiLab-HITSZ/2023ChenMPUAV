function [x,y1,y2,y3]=true_PS(Global_D,probID,t1,t2,varargin)
%%Calculte the approximate real common PS and PF forMPMOP
syms x;
g=1;
H=1000; % number of divisions along each objective.
if nargin==5
        t3=varargin{1};
end
y3=0;
D=Global_D;
switch (probID)
    case 'MPMOP1'
        a1=5*cos(0.5*pi*t1);
        a2=5*cos(0.5*pi*t2);
        x=linspace(1,4,20000);
        x1=find(abs(1./(1+exp(a1*(x-2.5)))-1./(1+exp(a2*(x-2.5))))<1e-4);
        temp=[];
        for i=1:(size(x1,2)-1)
            temp(i)=abs(x1(i)-x1(i+1))<10;
        end
        temp(size(x1,2))=0;
        x1=x1(:,~temp);
        x1=x(x1);
        xi=1/(1+exp(a2*(x1-2.5)));
        x=x1;
        for i=1:D-1
            x=[x xi];
        end     
        y1=MPMOP_Value('MPMOP1', x, t1);
        y2=MPMOP_Value('MPMOP1', x, t2);
        
    case 'MPMOP2'
        G1=sin(0.5*pi*t1);
        G2=sin(0.5*pi*t2);
        x=linspace(0,1,40000);
        x1=find(abs(G1*sin(4*pi*x)/(1+abs(G1))-G2*sin(4*pi*x)/(1+abs(G2)))<1e-4);
        temp=[];
        for i=1:(size(x1,2)-1)
            temp(i)=abs(x1(i)-x1(i+1))<10;
        end
        temp(size(x1,2))=0;
        x1=x1(:,~temp);
        x1=x(x1);
        xi=G2*sin(4*pi*x1)/(1+abs(G2));
        x=x1';
        for i=1:D-1
            x=[x xi'];
        end       
        y1=MPMOP_Value('MPMOP2', x, t1);
        y2=MPMOP_Value('MPMOP2', x, t2);
        
    case 'MPMOP3' 
        N1=1+floor(10*abs(sin(0.5*pi*t1)));
        N2=1+floor(10*abs(sin(0.5*pi*t2)));
        g=1;
        x1=linspace(0,1,10000);  
        PS=true(1,10000);
        for i=0:N1-1
            PS(find(x1>i/N1&x1<(2*i+1)/(2*N1)))=false;
        end
        for i=0:N2-1
            PS(find(x1>i/N2&x1<(2*i+1)/(2*N2)))=false;
        end
        x1=x1(PS);
        x=x1';
        for i=2:D
            xi=cos(4*t1+x(:,1)+x(:,end));
            x=[x xi];
        end   
        y1=MPMOP_Value('MPMOP3', x, t1);
        y2=MPMOP_Value('MPMOP3', x, t2);       
    
    case 'MPMOP4'
        [x1,x2]=meshgrid(linspace(0,1,20*H));
        G1=sin(0.5*pi*t1);
        G2=sin(0.5*pi*t2);
        temp=abs(sin(2*pi*(x1+x2))./(1+abs(G1))-sin(2*pi*(x1+x2))./(1+abs(G2)))<1e-4;
        for i=1:(size(x1,1)-1)
            for j=1:(size(x1,2)-1)
                if (temp(i,j)==1 &temp(i,j+1)==1)| (temp(i,j)==1 &temp(i+1,j)==1)
                    temp(i,j)=0;
                end
            if (temp(i,j+1)==1 &temp(i+1,j+1)==1)
                temp(i,j+1)=0;
            end
            end
        end
        for j=1:(size(x1,2)-1)
            if(temp(i+1,j)==1 &temp(i+1,j+1)==1)
                temp(i+1,j)=0;
            end
        end
        x1=x1(temp);
        x2=x2(temp);
        temp=[];%release the space
        xi=sin(2*pi*(x1+x2))./(1+abs(G2));
        x=[x1,x2];
        for i=1:D-2
            x=[x xi];
        end   
        y1=MPMOP_Value('MPMOP4', x, t1);
        y2=MPMOP_Value('MPMOP4', x, t2);
        
    case 'MPMOP5'
        [x1,x2]=meshgrid(linspace(0,1,H));
         G1=abs(sin(0.5*pi*t1));
         G2=abs(sin(0.5*pi*t2));
         temp=abs(0.5*G1*x1-0.5*G2*x1)<1e-4;
         for i=1:(size(x1,1)-1)
            for j=1:(size(x1,2)-1)
                if (temp(i,j)==1 & temp(i,j+1)==1)
                    temp(i,j)=0;
                end
            end
        end
        for j=1:(size(x1,2)-1)
            if(temp(i+1,j)==1 & temp(i+1,j+1)==1)
                temp(i+1,j)=0;
            end
        end
        x1=x1(temp);
        x2=x2(temp);
        xi=0.5*G2*x1;
        x=[x1,x2];
        for i=1:D-2
            x=[x xi];
        end   
        y1=MPMOP_Value('MPMOP5', x, t1);
        y2=MPMOP_Value('MPMOP5', x, t2);

        
    case 'MPMOP6'
        [x1,x2]=meshgrid(linspace(0,1,H));
        temp=abs(sin(t1*x1)-sin(t2*x1))<1e-4;
        x1=x1(temp);
        x2=x2(temp);
        xi=sin(t2*x1);
        x=[x1,x2];
        for i=1:D-2
            x=[x xi];
        end   
        y1=MPMOP_Value('MPMOP6', x, t1);
        y2=MPMOP_Value('MPMOP6', x, t2);

        
    case 'MPMOP7'
        t=[t1,t2,t3];
        a=5*cos(0.5*pi*t);
        x=linspace(1,4,40000);
        x1=find((abs(1./(1+exp(a(1)*(x-2.5)))-1./(1+exp(a(2)*(x-2.5))))<1e-4) ...,
            &(abs(1./(1+exp(a(1)*(x-2.5)))-1./(1+exp(a(3)*(x-2.5))))<1e-4));
        temp=[];
        for i=1:(size(x1,2)-1)
            temp(i)=abs(x1(i)-x1(i+1))<10;
        end
        temp(size(x1,2))=0;
        x1=x1(:,~temp);
        x1=x(x1);
        xi=1./(1+exp(a(2)*(x1-2.5)));
        x=x1;
        for i=1:D-1
            x=[x xi];
        end  
        y1=MPMOP_Value('MPMOP1', x, t(1));
        y2=MPMOP_Value('MPMOP1', x, t(2));
        y3=MPMOP_Value('MPMOP1', x, t(3));

        
    case 'MPMOP8'
        t=[t1,t2,t3];
        G=sin(0.5*pi*t);
        x=linspace(0,1,40000);
        x1=find((abs(G(1)*sin(4*pi*x)/(1+abs(G(1)))-G(2)*sin(4*pi*x)/(1+abs(G(2))))<1e-4) ...,
            &(abs(G(1)*sin(4*pi*x)/(1+abs(G(1)))-G(3)*sin(4*pi*x)/(1+abs(G(3))))<1e-4));
        temp=[];
        for i=1:(size(x1,2)-1)
            temp(i)=abs(x1(i)-x1(i+1))<10;
        end
        temp(size(x1,2))=0;
        x1=x1(:,~temp);
        x1=x(x1);
        xi=G(2)*sin(4*pi*x1)/(1+abs(G(2)));
        x=x1';
        for i=1:D-1
            x=[x xi'];
        end   
         y1=MPMOP_Value('MPMOP2', x, t(1));
        y2=MPMOP_Value('MPMOP2', x, t(2));
        y3=MPMOP_Value('MPMOP2', x, t(3));


    case 'MPMOP9'
        [x1,x2]=meshgrid(linspace(0,1,20*H));
        G1=sin(0.5*pi*t1);
        G2=sin(0.5*pi*t2);
        G3=sin(0.5*pi*t3);
        temp=(abs(sin(2*pi*(x1+x2))./(1+abs(G1))-sin(2*pi*(x1+x2))./(1+abs(G2)))<1e-4) ...,
            &(abs(sin(2*pi*(x1+x2))./(1+abs(G1))-sin(2*pi*(x1+x2))./(1+abs(G3)))<1e-4);
        for i=1:(size(x1,1)-1)
            for j=1:(size(x1,2)-1)
                if (temp(i,j)==1 &temp(i,j+1)==1)| (temp(i,j)==1 &temp(i+1,j)==1)
                    temp(i,j)=0;
                end
            if (temp(i,j+1)==1 &temp(i+1,j+1)==1)
                temp(i,j+1)=0;
            end
            end
        end
        for j=1:(size(x1,2)-1)
            if(temp(i+1,j)==1 &temp(i+1,j+1)==1)
                temp(i+1,j)=0;
            end
        end 
        x1=x1(temp);
        x2=x2(temp);
        temp=[];%release the space
        xi=sin(2*pi*(x1+x2))./(1+abs(G2));
        x=[x1,x2];
        for i=1:D-2
            x=[x xi];
        end   
         y1=MPMOP_Value('MPMOP4', x, t1);
        y2=MPMOP_Value('MPMOP4', x, t2);
        y3=MPMOP_Value('MPMOP4', x, t3);
 
    case 'MPMOP10'
        [x1,x2]=meshgrid(linspace(0,1,H));
         G1=abs(sin(0.5*pi*t1));
         G2=abs(sin(0.5*pi*t2));
         G3=abs(sin(0.5*pi*t3));
         temp=(abs(0.5*G1*x1-0.5*G2*x1)<1e-4)&(abs(0.5*G1*x1-0.5*G3*x1)<1e-4);
         for i=1:(size(x1,1)-1)
            for j=1:(size(x1,2)-1)
                if (temp(i,j)==1 &temp(i,j+1)==1)
                    temp(i,j)=0;
                end
            end
        end
        for j=1:(size(x1,2)-1)
            if(temp(i+1,j)==1 &temp(i+1,j+1)==1)
                temp(i+1,j)=0;
            end
        end  
        x1=x1(temp);
        x2=x2(temp);
        xi=0.5*G2*x1;
        x=[x1,x2];
        for i=1:D-2
            x=[x xi];
        end   
         y1=MPMOP_Value('MPMOP5', x, t1);
        y2=MPMOP_Value('MPMOP5', x, t2);
        y3=MPMOP_Value('MPMOP5', x, t3);
 
    case 'MPMOP11'
        [x1,x2]=meshgrid(linspace(0,1,H));
        temp=(abs(sin(t1*x1)-sin(t2*x1))<1e-4)&(abs(sin(t1*x1)-sin(t3*x1))<1e-4);
        x1=x1(temp);
        x2=x2(temp);
        xi=sin(t2*x1);
        x=[x1,x2];
        for i=1:D-2
            x=[x xi];
        end   
         y1=MPMOP_Value('MPMOP6', x, t1);
        y2=MPMOP_Value('MPMOP6', x, t2);
        y3=MPMOP_Value('MPMOP6', x, t3);

       
   otherwise
        disp('no such test problem.')     
end
end


