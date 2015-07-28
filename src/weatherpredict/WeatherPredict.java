
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package weatherpredict;

/**
 *
 * @author ngo
 */
public class WeatherPredict {

    /**
     * @param args the command line arguments
     */
    private static final int X=50;
    private static final int Y=40;
    private static final int K=5;
    private static final int R=287;
    private static final int CP=1005;
    private static final int M=1;
    private static final int N=50;
    private static final int[] PN= {1000, 900, 650, 400, 150};
    private static final int P0=1000;
    private static final float delta=1200;
    private static final float dex=100000;
    private static final float F=(float) (14.584*Math.pow(10, -5));
    
    static float[][][] Ut= new float[X][Y][K];
    static float[][][] Ut1= new float[X][Y][K];
    static float[][][] Ut2= new float[X][Y][K];
    
    static float[][][] Vt= new float[X][Y][K];
    static float[][][] Vt1= new float[X][Y][K];
    static float[][][] Vt2= new float[X][Y][K];
    
    static float[][][] TE_t= new float[X][Y][K];
    static float[][][] TE_t1= new float[X][Y][K];
    static float[][][] TE_t2= new float[X][Y][K];
    static float[][][] Q_t= new float[X][Y][K];
    static float[][][] Q_t1= new float[X][Y][K];
    static float[][][] Q_t2= new float[X][Y][K];
    static float[][] PS= new float[X][Y];
    static float[][] PS_t1= new float[X][Y];
    static float[][] PS_t2= new float[X][Y];
    static float[][][] Wt= new float[X][Y][K];
    static float[][][] Ws= new float[X][Y][1];
    static float[][][] PHI= new float[X][Y][K];
    
    public void init(){
        
    }
    public static void main(String[] args) {
        // TODO code application logic here
        float M1= M/(4*deltata);
        float RC= R/CP;
        float P1, P2, P3, P4; 
        float W9;
        float R0, TV,V_TB;
        for(int n=0; n<N; n++){
            for(int k=0; k<K; k++){
                //giai doan du doan
                for (int i=1;i<X-1; i++){
                    for(int j=1; j<Y-1; j++){
                        //2.1
                        Ut2[i][j][k]=-M*(PHI[i+1][j][k]-PHI[i][j][k])/dex
                            -M1*((Ut1[i][j][k]+Ut1[i+1][j][k])*(Ut1[i+1][j][k]-Ut1[i][j][k])+(Ut1[i][j][k]+Ut1[i-1][j][k])*(Ut1[i][j][k]-Ut1[i-1][j][k]))
                            +F*(Vt1[i][j][k]+Vt1[i][j-1][k]+Vt1[i+1][j][k]+Vt1[i+1][j-1][k])/4
                            -M1*((Vt1[i][j][k]+Vt1[i+1][j][k])*(Ut1[i][j+1][k]-Ut1[i][j][k])+(Vt1[i][j-1][k]+Vt1[i+1][j-1][k])*(Ut1[i][j-i][k]-Ut1[i][j][k]));
                        //2.2
                        Vt2[i][j][k]= -M*(PHI[i][j+1][k]-PHI[i][j][k])/deltata
                            - M1*((Ut1[i][j][k]+Ut1[i][j+1][k])*(Vt1[i+1][j][k]-Vt1[i][j][k])+(Ut1[i-1][j][k]+Ut1[i-1][j+1][k])*(Vt1[i][j][k]-Vt1[i][j][k]))
                            - M1*((Vt1[i][j][k]+Vt1[i][j+1][k])*(Vt1[i][j+1][k]-Vt1[i][j][k])+(Vt1[i][j][k]+Vt1[i][j-1][k])*(Vt1[i][j][k]-Vt1[i][j-1][k]))
                            -F*(Ut1[i][j][k]+Ut1[i-1][j][k]+Ut1[i][j+1][k]+Vt1[i-1][j+1][k])/4;
                        //2.3
                        TE_t2[i][j][k]=-2*M1*(Ut1[i+1][j][k]*(TE_t1[i+1][j][k]-TE_t1[i][j][k])
                                -Ut1[i][j][k]*(TE_t1[i][j][k]-TE_t1[i-1][j][k])
                                +Vt1[i][j+1][k]*(TE_t1[i][j+1][k]-TE_t1[i][j][k])
                                -Vt1[i][j][k]*(TE_t1[i][j][k]-TE_t1[i][j-1][k]));
                        //2.4
                        Q_t2[i][j][k]=-2*M1*(Ut1[i][j][k]*(Q_t1[i+1][j][k]-Q_t1[i][j][k])
                                -Ut1[i][j][k]*(Q_t1[i][j][k]-Q_t1[i-1][j][k])
                                +V[n-1][i][j+1][k]*(Q_t1[i][j+1][k]-Q_t1[i][j][k])
                                -Vt1[i][j][k]*(Q_t1[i][j][k]-Q_t1[i][j-1][k]));
                        if(k==0){
                            P1= (PN[k+1]+PN[k])/2;                            
                            P2= (PN[k+2]+PN[k+1])/2;                                                        
                            //**1
                            Ut2[i][j][k]=Ut2[i][j][k]-(Wt[i+1][j][k]+Wt[i][j][k])*(Ut1[i][j][k+1]-Ut1[i][j][k])/(2*(PN[k+1]-PN[k]));
                            
                            //**2
                            Vt2[i][j][k]=Vt2[i][j][k]-(Wt[i][j+1][k]+Wt[i][j][k])*(Vt1[i][j][k+1]-Vt1[i][j][k])/(2*(PN[k+1]-PN[k]));
                            
                            //**3
                            TE_t2[i][j][k]=TE_t2[i][j][k]-Wt[i][j][k]*(TE_t1[i][j][k+1]-TE_t1[i][j][k])/(PN[k+1]-PN[k]);

                            //**4
                            Q_t2[i][j][k]= Q_t2[i][j][k]-Wt[i][j][k]*(Q_t1[i][j][k+1]-Q_t1[i][j][k])/(PN[k+1]-PN[k]);

                            //**5
                            PS_t2[i][j]=-2*M1*(Ut1[i+1][j][k]*(PS_t1[i+1][j]-PS_t1[i][j])
                                    -Ut1[i][j][k]*(PS_t1[i][j]-PS_t1[i-1][j])
                                    +Vt1[i][j+1][k]*(PS_t1[i][j+1]-PS_t1[i][j]));

                            W9= (Wt[i][j][k+2]*P1+Wt[i][j][k]*P2)/(P1+P2);
                            PS_t2[i][j]=PS_t2[i][j]+W9-(M*((Ut1[i+1][j][k+1]-Ut1[i][j][k])
                                    +(Vt1[i][j+1][k+1]-Vt1[i][j][k+1])+(Ut1[i+1][j][k]-Ut1[i][j][k])
                                    +(Vt1[i][j+1][k]-Vt1[i][j][k]))/(2*delta))*(PS_t1[i][j]-900);
                        }else if(k==K-1){
                            //***
                            P1= (PN[k]-PN[k-1])/2;                            
                            P2=25;                            
                            //***1
                            Ut2[i][j][k]= Ut2[i][j][k]
                                    -(Wt[i+1][j][k-1]+Wt[i][j][k-1])*(PN[k]-P2)/((2*(P1- P2)))*((Ut1[i][j][k]-Ut1[i][j][k-1])/(PN[k]-PN[k-1]));

                            //***2
                            Vt2[i][j][k]= Vt2[i][j][k]
                                    -((Wt[i][j+1][k]+Wt[i][j][k-1])*(PN[k]-PN[k-1])/(2*(P1-P2)))*((Vt1[i][j][k-1]-Vt1[i][j][k-1])/(PN[k]-PN[k-1]));

                            //***3
                            TE_t2[i][j][k]=0;
                            Q_t2[i][j][k]=0;

                        }else{
                        P1= (PN[k-1]+PN[k])/2;
                        P2= (PN[k]+PN[k+1])/2;
                        //2.5
                        Ut2[i][j][k]= Ut2[i][j][k]-(((Wt[i+1][j][k]+Wt[i][j][k])*(PN[k]-P2)
                                +(Wt[i+1][j][k+1]+Wt[i][j][k+1])*(P1-PN[k]))/(2*(P1-P2))
                                *((Ut1[i][j][k-1]-Ut1[i][j][k])*(PN[k]-PN[k+1])/(PN[k-1]-PN[k]))
                                +((Ut1[i][j][k]-Ut1[i][j][k+1])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1])))/(PN[k-1]-PN[k+1]);
                        //2.6
                        Vt2[i][j][k]= Vt2[i][j][k]-(((Wt[i][j+1][k]+Wt[i][j][k])*(PN[k]-P2)
                                +(Wt[i][j+1][k+1]+Wt[i][j][k+1])*(P1-PN[k]))/(2*(P1-P2))
                                *((Vt1[i][j][k-1]-Vt1[i][j][k+1])*(PN[k]-PN[k+1])/(PN[k-1]-PN[k]))
                                +((Vt1[i][j][k]-Vt1[i][j][k+1])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1])))/(PN[k-1]-PN[k+1]);
                        //2.7
                        TE_t2[i][j][k]=TE_t2[i][j][k]-Wt[i][j][k]*((TE_t1[i][j][k-1]-TE_t1[i][j][k])
                                *(PN[k]-PN[k+1])/(PN[k-1]-PN[k])+(TE_t1[i][j][k]-TE_[n][i][j][k+1])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1]))/(PN[k-1]-PN[k+1]);
                        //2.8
                        Q_t2[i][j][k]= Q_t2[i][j][k]-Wt[i][j][k]*((Q_t1[i][j][k-1]-Q_t1[i][j][k])*
                                (PN[k]-PN[k+1])/(PN[k-1]-PN[k])+(Q_t1[i][j][k]-Q_t1[i][j][k+1])*
                                (PN[k-1]-PN[k])/(PN[k]-PN[k+1]))/(PN[k-1]-PN[k+1]);

                        }
                    }
                }
                //giai doan chuan doan
                for (int i=1; i<X-1; i++){
                    for(int j=1; j<Y-1; j++){
                        R0= (float) (TE_t1[i][j][k]*Math.pow(PN[k]/P0, RC));
                        TV=(float) (TE_t1[i][j][k]*(1+0.61*Q_t1[i][j][k]));
                        if(k==0){
                            P1=(PN[k]-PN[k+1])/2;
                            P2=(PN[k+1]-PN[k+2])/2;
                            P3=(PN[k+1]+PN[k+2])/2;
                            P4=(PN[1]+PN[2])/2;
                            //**
                            V_TB=-M*(((Ut1[i+1][j][k]-Ut1[i][j][k])
                                 +(Vt1[i][j+1][k]-Vt1[i][j][k]))*P2
                                 +((Ut1[i+1][j][k+1]-Ut1[i][j][k+1])
                                 +(Vt1[i][j+1][k]-Vt1[i][j][k+1]))*P1)/(delta*(P1+P2));
                            float VTS= M*((Ut1[i+1][j][k]-Ut1[i][j][k])+(Vt1[i][j+1][k]-Vt1[i][j][k]))/delta;
                            Ws[i][j][0]=Ws[i][j][0]+VTS*(P4- PN[k]);
                            //**1
                            Wt[i][j][k]=(-V_TB*(PN[k]-P3)+(Ws[i][j][1]*((PN[k+1]-P3)/(PN[k]-PN[k+1])))
                                    -(Wt[i][j][k+1]*((PN[k]-PN[k+1])/(PN[k+1]-P3))))/((PN[k+1]-P3)/(PN[k]-PN[k+1])-(PN[k]-PN[k+1])/(PN[k+1]-P3));
                            //**2
                            PHI[i][j][k]=(float) (PHI[i][j][k+1]+R/R0*TV*Math.pow((PN[k]/P0), RC)*(PN[k+1]-PN[k]));
                        } else if (k==K-1){
                            //***1
                            PHI[i][j][k]= (float) (PHI[i][j][k-1]+R/R0*TV*Math.pow(PN[k]/P0, RC)*(PN[k]-PN[k-1]));
                            Wt[i][j][k]=0;
                        }else{
                            PHI[i][j][k]=(float) ((R/R0*TV*Math.pow(PN[k]/P0, RC)*(PN[k-1]-PN[k+1])+PHI[i][j][k+1]*((PN[k]-PN[k+1])/(PN[k-1]-PN[k]))
                                    -PHI[i][j][k+1]*((PN[k-1]-PN[k])/(PN[k]-PN[k-1])))/((PN[k]-PN[k+1])/(PN[k-1]-PN[k])-((PN[k-1]-PN[k])/(PN[k]-PN[k+1]))));
                            V_TB=-M*((Ut1[i+1][j][k]-Ut1[i][j][k])
                                    +(Vt1[i][j+1][k]-Vt1[i][j][k])
                                    +(Ut1[i+1][j][k+1]-Ut1[i][j][k+1])
                                    +(Vt1[i][j+1][k+1]-Vt1[i][j][k+1]))/(2*delta);
                            Wt[i][j][k]=(-V_TB*(PN[k-1]-PN[k+1])
                                    +Wt[i][j][k-1]*((PN[k]-PN[k+1])/(PN[k-1]-PN[k]))
                                    -Wt[i][j][k+1]*((PN[k-1]-PN[k])/(PN[k]-PN[k+1])))
                                    /((PN[k]-PN[k+1])/(PN[k-1]-PN[k])-(PN[k-1]-PN[k])/(PN[k]-PN[k+1]));
                        }
                        
                    }
                }
            }
        }
    }
}
