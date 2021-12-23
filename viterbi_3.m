%This is a soft decision based viterbi decoder for constraint length 3.
%We are encoding 300 data bits using the convolutional encoder block
%available in matlab.
%The code is 1/2 rate encoded with generator polynomials (7,5). 
%Hence, after decoding, we will get half of 600 bits, i.e. 300 bits, which represents the original information. */

input=[];

data=randi([0 1],18,1);%generating random data bits
d2=data;


trellis=poly2trellis(3,[7 5]);
codedata=convenc(data,trellis);;%convolutional encoding
cd=codedata';

k=1;

for i=1:1:18
      for j=1:1:2
             bits(i,j)=cd(1,k);% performing grouping of the encoded bits into a 300 by 2 matrix for computation.
             k=k+1;
      end
end
    
% Each convolution encoder has a finite number of possible states depending on constraint length.
%number of possible states=2^(constraint length-1)
input=0:3; %matrix storing all possible states
data=decimalToBinaryVector(input);

%calculating G1 G2 outputs
outputstate0g1=[];
outputstate0g2=[];
for i=1:4
   outputstate0g1(i)=xor(data(i,1),data(i,2));
   outputstate0g2(i)=data(i,2);
end
%corresponding to each set of state transition there is an output. 
%out0 stores output corresponding to transition from n0 to data and similarly out1 for transition from n1 to data.


os0g1=(outputstate0g1)';
os0g2=(outputstate0g2)';
out0=horzcat(os0g1,os0g2);
out1=~(out0); 


%initialising metrics table
metrics=ones(4,18);
metrics=-1*metrics;

predecessor=ones(4,4);
predecessor=predecessor*(-1);
predecessor(1,1)=0; %intialising in all zero state
flag=-1;
%Each current state of the convolution encoder will have two possible previous states, stored in n0 and n1. 
%n0 represents that state with LSB 0, while n1 with LSB=1.

ns=data;
n0=circshift(ns,[0,1]); 
n0(:,2)=0;
n1=circshift(ns,[0,1]);
n1(:,2)=1;

%computing predecessor table till constraint length(=3)
for i=2:1:3
    for j=1:1:4
        for k=1:1:4
            if (n0(j,:)==data(k,:))%capturing the index of n0 using the table of data
                n0d=k    ;     %n0d stores captured index                   
            end
        end
        
        for m=1:1:4
            if(n1(j,:)==data(k,:))%capturing the index of n1 using the table of data
                n1d=m;%n1d stores capured index
            end
        end
       %till the depth equal to the constraint length, each current state only has one predecessor state.
 %If it is a valid  predecessor, then it will have a positive value, else -1.          
                
        if (~(predecessor(n0d,(i-1))==flag) | ~(predecessor(n1d,(i-1))==flag))
          %  d=3
            if(predecessor(n0d,(i-1))~=flag)
            predecessor(j,i)=n0d-1;
            else
                predecessor(j,i)=n1d-1;
            end
        end
    end
    
end

%computing metrics table till constraint length (=3)
metrics(1,1)=0;
for i=2:1:3
    for j=1:1:4
       x=predecessor(j,i);
       if(~(x==-1))
           msb=data(j,1);
           
           if(msb==0)%/this implies an input of 0 caused this state transition, hence refer to out0
               t0(1,1)=out0((x+1),1);
               t0(1,2)=out0((x+1),2);
           end 
           
           if(msb==1)%this implies an input of 1 caused this state transition, hence refer to out1
               t0(1,1)=out1((x+1),1);
               t0(1,2)=out1((x+1),2);
           end 
           %accumulated metric is distance between the i/p and o/p bits + previous stage metric
           metrics(j,i)=xor(t0(1,1),bits((i-1),1))+xor(t0(1,2),bits((i-1),2))+metrics((x+1),(i-1));
           
       end
    end
end 

%computing metrics and predecessor table for the rest of the time
%instants(t>=7)
%each state has 2 possible predecessor states, of which the one with larger accumulated metric is to be eliminated.
%the two possible metrics are stored in a and b.

for k=3:1:17
for i=1:1:4
    x=ns(i,1);%storing the digit which must have been input to produce the nextstate
        for j=1:1:4
            if (n0(i,:)==data(j,:))
                 if(x==0)%if '0' was input, corresponding output has been stored
                     y0=out0(j,:);
                 end
                 if(x==1)%if '1' was input, corresponding output has been stored
                     y0=out1(j,:);
                 end
                 %first possible accumulated metric
                
                 a=sqrt((bits(k,1)-y0(1,1))^2+(bits(k,2)-y0(1,2))^2)+metrics(j,k);
                 inx1=j;
            end
            
        end    
         
        
        for j=1:1:4 
            if (n1(i,:)==data(j,:))
                 if(x==0)%if '0' was input, corresponding output has been stored
                     y1=out0(j,:);
                 end
                 if(x==1)%if '1' was input, corresponding output has been stored
                     y1=out1(j,:);
                 end
                  %second possible accumulated metric
                  b=sqrt((bits(k,1)-y1(1,1))^2+(bits(k,2)-y1(1,2))^2)+metrics(j,k);
                  inx2=j;%capturing the predecessor state
            end           
        end
       
        metrics(i,k+1)= min(a,b);
         % comparing the two to eliminate the one with larger metric.
         if(min(a,b)==a)
            predecessor(i,k+1)=inx1-1;
        else
            predecessor(i,k+1)=inx2-1;
        end
        
end       
end 

%TRACEBACK MECHANISM
%After the construction of predecessor state stable,
%tracebacking through the trellis occurs and 'states' stores the respective state for the path to be taken for decoding.
%making the state row vector

val=predecessor(1,18);
for j=18:-1:2
    states(j-1)=val;
    val=predecessor((val+1),(j-1));
end

%corresponding to each state transition, we have an input bit.
%This is captured and stored in decoded.


for i=1:1:4
    
    decdata(i)=data(i,1)*2+data(i,2);
end

for i=2:1:17
    for j=1:1:4
        if states(i)==decdata(j)
            decoded(i-1)=data(j,1);
        end
    end
end

