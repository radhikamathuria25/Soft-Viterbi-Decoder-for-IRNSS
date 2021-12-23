/* This is a soft decision based viterbi decoder for constraint length 3. 
We have taken some arbitrary encoded bits which will be feeded in the decoder. The code is 1/2 rate encoded with generator polynomials (7,5). 
Hence, after decoding, we will get half of 36 bits, i.e. 18 bits, which represents the original information. */


#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>


int main()
{
    int bits[18][2],binary[4][2],n0[4][2],n1[4][2],data[4][2],out0[4][2],out1[4][2],pre[4][18],ns[4][2];


    int cd[36]={0,0 ,1 ,1 ,1,1,0,0,0,1,1,0,0,1,1,1,1,1,1,0,0,0,0,0,1,1,0,0,1,1,1,0,1,1,0,0};//taking some arbitrary convolutionally encoded bits

    int k=0;
    for (int i=0;i<18;i++)
    {       for (int j=0;j<2;j++)
        {    bits[i][j]=cd[k];// performing grouping of the encoded bits into a 18 by 2 matrix for computation.
            k++;
                   
        }
        
    }
    
    // Each convolution encoder has a finite number of possible states. All possible states are produced and stored in data.
    
    for(int a=0;a<4;a++)
    {   int b[2],b2[2];
       for (int i = 1; i >= 0; i--) 
    { 
        k = a>>i; 
        if (k & 1) 
          b[i]=1;
        else
            b[i]=0;
    } 
    for (int i=1;i>=0;i--)
    {   b2[1-i]=b[i]; 
       data[a][1-i]=b2[1-i];
    }
   
    }
    


    //Each current state of the convolution encoder will have two possible previous states, stored in n0 and n1. n0 represents that state with LSB 0, while n1 with LSB=1.
    for(int i=0;i<4;i++)
    {
            for (int j=0;j<1;++j)
            {
                n0[i][j]=data[i][j+1];
                n0[i][1]=0;
                n1[i][j]=data[i][j+1];
                n1[i][1]=1;
                ns[i][j]=data[i][j];

            }
    }
    
    //corresponding to each set of state transition there is an output. out0 stores output corresponding to transition from n0 to data and similarly out1 for transition from n1 to data.

   for (int i=0;i<4;i++)
   {
       out0[i][0]=data[i][0]^data[i][1];
       out1[i][0]=1^out0[i][0];
       out0[i][1]=data[i][1];
       out1[i][1]=1^out0[i][1];
   }
    
    //COMPUTING PREDECESSOR TILL CONSTRAINT LENGTH

    for(int i=0;i<4;i++)
    {
        for(int j=0;j<3;j++)
        {
            pre[i][j]=-1;
        }
    }
    pre[0][0]=0;
    int flag=-1,cnt;
    int n0d,n1d;

    for (int i=1;i<3;i++)
    {
        for (int j=0;j<4;j++)
            {
                for(int r=0;r<4;r++)
                { cnt=0;
                    for(int c=0;c<2;c++)
                    {
                        if(n0[j][c]!=data[r][c]) //capturing the index of n0 in the table of data 
                        {
                            cnt++;
                            break;
                        } 
                    }
                    if(cnt==0)
                    n0d=r;//n0d stores captured index 

                }

                for(int r=0;r<4;r++)
                {   cnt=0;
                    for(int c=0;c<2;c++)
                    {
                        if(n1[j][c]!=data[r][c])//capturing the index of n1 in the data table
                        { 
                            cnt++;
                            break;
                        }
                    }
                    if(cnt==0)
                    n1d=r;//n1d stores captured index
                }

                if((pre[n1d][i-1]!=flag) || (pre[n0d][i-1]!=flag))//till the depth equal to the constraint length, each current state only has one predecessor state. If it is a valid  predecessor, then it will have a positive value, else -1.
                {
                    if(pre[n0d][i-1]!=flag)
                    pre[j][i]=n0d;
                    else
                    pre[j][i]=n1d;
                }
            }
    }
        
//computing metrics table till constraint length

int metrics[4][18];

//initialising entire metrics to -1
for (int i=0;i<4;i++)
{
    for(int j=0;j<18;j++)
    {
        metrics[i][j]=-1;
    }
}

metrics[0][0]=0;
int x,msb,t0[1][2];

for(int i=1;i<=2;i++)
{
    for(int j=0;j<4;j++)
    {
       x=pre[j][i];

       if(x!=-1)
       {
           msb=data[j][0];

           if(msb==0)//this implies an input of 0 caused this state transition, hence refer to out0
           {
               t0[0][0]=out0[x][0];
               t0[0][1]=out0[x][1];
           }

            if(msb==1)//this implies an input of 1 caused this state transition, hence refer to out1
           {
               t0[0][0]=out1[x][0];
               t0[0][1]=out1[x][1];
           }

          

        metrics[j][i]=(t0[0][0]^bits[(i-1)][0] )+ ( t0[0][1]^bits[(i-1)][1] )+ metrics[x][(i-1)];//accumulated metric is distance between the i/p and o/p bits + previous stage metric
        
       } 
    }
}


//computing rest of the metrics and predecessor states. capturing the index has been done as above.
//each state has 2 possible predecessor states, of which the one with larger accumulated metric is to be eliminated. the two possible metrics are stored in a and b.
int y,inx1,inx2,y0[1][2],y1[1][2],a,b;
for(int m=2; m<18; m++)
{   
    for(int i=0;i<4;i++)
    {
        y=data[i][0];
       for (int j=0;j<4;j++)
        {
            cnt=0;
            for (int c=0;c<2;c++)
            {
                if(n0[i][c]!=data[j][c])
                {
                     cnt++;
                    break;
                }
            }

            if(cnt==0)
            {
                if(y==0)//this implies an input of 0 caused this state transition, hence refer to out0
                {
                    y0[0][0]=out0[j][0];
                    y0[0][1]=out0[j][1];
                }

                if(y==1)//this implies an input of 1 caused this state transition, hence refer to out1       
                {
                    y0[0][0]=out1[j][0];
                    y0[0][1]=out1[j][1];
                }
                a=(pow((bits[m][0]-y0[0][0]),2))+(pow((bits[m][1]-y0[0][1]),2))+metrics[j][m];//first possible accumulated metric
                inx1=j;
            }
        }

        for (int j=0;j<4;j++)
        { 
             cnt=0;
            for (int c=0;c<2;c++)
            {
                if(n1[i][c]!=data[j][c])
                {
                    cnt++;
                    break;
                }
            }

            if(cnt==0)
            {
                if(y==0)
                {
                    y1[0][0]=out0[j][0];
                    y1[0][1]=out0[j][1];
                }

                if(y==1)
                {
                    y1[0][0]=out1[j][0];
                    y1[0][1]=out1[j][1];
                }

                b=(pow((bits[m][0]-y1[0][0]),2))+(pow((bits[m][1]-y1[0][1]),2))+metrics[j][m];//second possible accumulated metric
                inx2=j;//capturing the predecessor state
            }

        }
        int min;
        if(a>b)//comparing the two to eliminate the one with larger metric.
        min=b;
        else
        min=a;
        metrics[i][m+1]=min;
        if(min==a)
        pre[i][m+1]=inx1;
        else
        pre[i][m+1]=inx2;
    }

}

//TRACEBACK MECHANISM
//After the construction of predecessor state stable, tracebacking through the trellis occurs and 'states' stores the respective state for the path to be taken for decoding.
int val=pre[0][17];
int states[17],decoded[17];
for (int j=17;j>=1;j--)
{
    states[j-1]=val;
    val=pre[val][j-1];

}
//corresponding to each state transition, we have an input bit. This is captured and stored in decoded.
for(int i=1;i<18;i++)
{
    for(int j=0;j<4;j++)
    {
        if(states[i]==j)
        decoded[i-1]=data[j][0];
    }

}

//printing the final decoded 16 bits
    for(int j=0;j<16;j++)
    {
        printf("%d",decoded[j]);
    }
    printf("\n");

    return 0;
}

