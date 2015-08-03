
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package weatherpredict;

import java.util.Scanner;

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
    private static final float det=1200;
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
    
    static float P1, P2, P3, P4, P5; 
    static float U1, U2, V1, V2;
    static float W1, W2, W9;
    static float R0, TV,V_TB;
    static float M1= M/(4*dex);
    static float RC= R/CP;
        
    static float det1=2*det;
    static Scanner scanner;
    static boolean isFinished;
    public static void init(){
        System.out.println("Khoi tao so lieu");
        System.out.println("1. Nhap bang tay");
        System.out.println("2. Doc tu file");
        int option= scanner.nextInt();
        if(option==1) initManually();
        else if(option==2) initByFile();
        isFinished=false;
    }
    //nhap du lieu bang tay
    public static void initManually(){
        
    }
    //nhap du lieu tu file
    public static void initByFile(){
        
    }
    public static void DuDoan(){
        for(int k=0; k<K; k++){                
            for (int i=1;i<X-1; i++){
                for(int j=1; j<Y-1; j++){
                    //2.1
                    Ut2[i][j][k]=-M*(PHI[i+1][j][k]-PHI[i][j][k])/dex
                        -M1*((Ut1[i][j][k]+Ut1[i+1][j][k])*(Ut1[i+1][j][k]-Ut1[i][j][k])+(Ut1[i][j][k]+Ut1[i-1][j][k])*(Ut1[i][j][k]-Ut1[i-1][j][k]))
                        +F*(Vt1[i][j][k]+Vt1[i][j-1][k]+Vt1[i+1][j][k]+Vt1[i+1][j-1][k])/4
                        -M1*((Vt1[i][j][k]+Vt1[i+1][j][k])*(Ut1[i][j+1][k]-Ut1[i][j][k])+(Vt1[i][j-1][k]+Vt1[i+1][j-1][k])*(Ut1[i][j-i][k]-Ut1[i][j][k]));
                    //2.2
                    Vt2[i][j][k]= -M*(PHI[i][j+1][k]-PHI[i][j][k])/dex
                        - M1*((Ut1[i][j][k]+Ut1[i][j+1][k])*(Vt1[i+1][j][k]-Vt1[i][j][k])+(Ut1[i-1][j][k]+Ut1[i-1][j+1][k])*(Vt1[i][j][k]-Vt1[i][j][k]))
                        - M1*((Vt1[i][j][k]+Vt1[i][j+1][k])*(Vt1[i][j+1][k]-Vt1[i][j][k])+(Vt1[i][j][k]+Vt1[i][j-1][k])*(Vt1[i][j][k]-Vt1[i][j-1][k]))
                        -F*(Ut1[i][j][k]+Ut1[i-1][j][k]+Ut1[i][j+1][k]+Vt1[i-1][j+1][k])/4;
                    U1= (Ut1[i+1][j][k]+Ut1[i+1][j][k+1])/2;
                    U2= (Ut1[i][j][k]+Ut1[i][j][k+1])/2;
                    V1= (Vt1[i][j+1][k]+Vt1[i][j+1][k+1])/2;
                    V2=(Vt1[i][j][k]+Vt1[i][j][k+1])/2;
                    //2.3
                    TE_t2[i][j][k]=-2*M1*(U1*(TE_t1[i+1][j][k]-TE_t1[i][j][k])
                            -U2*(TE_t1[i][j][k]-TE_t1[i-1][j][k])
                            +V1*(TE_t1[i][j+1][k]-TE_t1[i][j][k])
                            -V2*(TE_t1[i][j][k]-TE_t1[i][j-1][k]));
                    //2.4
                    Q_t2[i][j][k]=-2*M1*(U1*(Q_t1[i+1][j][k]-Q_t1[i][j][k])
                            -U2*(Q_t1[i][j][k]-Q_t1[i-1][j][k])
                            +V1*(Q_t1[i][j+1][k]-Q_t1[i][j][k])
                            -V2*(Q_t1[i][j][k]-Q_t1[i][j-1][k]));
                    if(k==0){
                        P1= (PN[k+1]+PN[k])/2;                            
                        P2= (PN[k+2]+PN[k+1])/2;    
                        W9= (Wt[i][j][k+2]*P1+Wt[i][j][k]*P2)/(P1+P2);
                        W1= ((Wt[i+1][j][k]+Wt[i][j][k])+(Ws[i][j][0]+Ws[i+1][j][0]))/4;
                        W2=((Wt[i][j+1][k]+Wt[i][j][k])+(Ws[i][j+1][0]+Ws[i][j][0]))/4;
                        //**1
                        Ut2[i][j][k]=Ut2[i][j][k]-W1*(Ut1[i][j][k+1]-Ut1[i][j][k])/(PN[k+1]-PN[k]);

                        //**2
                        Vt2[i][j][k]=Vt2[i][j][k]-W2*(Vt1[i][j][k+1]-Vt1[i][j][k])/(PN[k+1]-PN[k]);

                        //**3
                        TE_t2[i][j][k]=TE_t2[i][j][k]-Wt[i][j][k]*(TE_t1[i][j][k+1]-TE_t1[i][j][k])/(PN[k+1]-PN[k]);

                        //**4
                        Q_t2[i][j][k]= Q_t2[i][j][k]-Wt[i][j][k]*(Q_t1[i][j][k+1]-Q_t1[i][j][k])/(PN[k+1]-PN[k]);

                        //**5
                        PS_t2[i][j]=-2*M1*(Ut1[i+1][j][k]*(PS_t1[i+1][j]-PS_t1[i][j])
                                -Ut1[i][j][k]*(PS_t1[i][j]-PS_t1[i-1][j])
                                +Vt1[i][j+1][k]*(PS_t1[i][j+1]-PS_t1[i][j]));


                        PS_t2[i][j]=PS_t2[i][j]+W9-(M*((Ut1[i+1][j][k+1]-Ut1[i][j][k])
                                +(Vt1[i][j+1][k+1]-Vt1[i][j][k+1])+(Ut1[i+1][j][k]-Ut1[i][j][k])
                                +(Vt1[i][j+1][k]-Vt1[i][j][k]))*2*M1)*(PS_t1[i][j]-900);
                    }else if(k==K-1){
                        //***
                        P1= (PN[k]+PN[k-1])/2;                            
                        P2=25;  
                        P5=(PN[k]-P2)/(2*(P1-P2));
                        //***1
                        Ut2[i][j][k]= Ut2[i][j][k]
                                -((Wt[i+1][j][k-1]+Wt[i][j][k-1])*P5*(Ut1[i][j][k]-Ut1[i][j][k-1]))/(PN[k]-PN[k-1]);

                        //***2
                        Vt2[i][j][k]= Vt2[i][j][k]
                                -((Wt[i][j+1][k]+Wt[i][j][k-1])*P5*(Vt1[i][j][k-1]-Vt1[i][j][k-1]))/(PN[k]-PN[k-1]);

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
                            *(PN[k]-PN[k+1])/(PN[k-1]-PN[k])+(TE_t1[i][j][k]-TE_t[i][j][k+1])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1]))/(PN[k-1]-PN[k+1]);
                    //2.8
                    Q_t2[i][j][k]= Q_t2[i][j][k]-Wt[i][j][k]*((Q_t1[i][j][k-1]-Q_t1[i][j][k])*
                            (PN[k]-PN[k+1])/(PN[k-1]-PN[k])+(Q_t1[i][j][k]-Q_t1[i][j][k+1])*
                            (PN[k-1]-PN[k])/(PN[k]-PN[k+1]))/(PN[k-1]-PN[k+1]);

                    }

                }
            }
        }
        //xu li ket qua du doan
        for (int k = 0; k < K; k++) {
            for (int i = 0; i < X; i++) {
                for (int j = 0; j < Y; j++) {
                    Ut2[i][j][k]= Ut[i][j][k]+det1*Ut2[i][j][k];
                    Vt2[i][j][k]= Vt[i][j][k]+det1*Vt2[i][j][k];
                    TE_t2[i][j][k]= TE_t[i][j][k]+det1*TE_t2[i][j][k];
                    Q_t2[i][j][k]= Q_t[i][j][k]+det1*Q_t2[i][j][k];                        
                    Ut[i][j][k]=Ut1[i][j][k]; Ut1[i][j][k]=Ut2[i][j][k];
                    Vt[i][j][k]=Vt1[i][j][k]; Vt1[i][j][k]=Vt2[i][j][k];
                    TE_t[i][j][k]=TE_t1[i][j][k]; TE_t1[i][j][k]=TE_t2[i][j][k];
                    Q_t[i][j][k]=Q_t1[i][j][k]; Q_t1[i][j][k]=Q_t2[i][j][k];
                }
            }
        }
        for (int i = 1; i < X-1; i++) {
            for (int j = 1; j < Y-1; j++) {
                PS_t2[i][j]= PS[i][j]+det1*PS_t2[i][j];
                PS[i][j]=PS_t1[i][j];
                PS_t1[i][j]=PS_t2[i][j];
            }                
        }
    }
    
    public static void ChuanDoan(){
        for(int k=0; k<K; k++){
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
                            +(Vt1[i][j+1][k]-Vt1[i][j][k+1]))*P1)/(dex*(P1+P2));
                        float VTS= M*((Ut1[i+1][j][k]-Ut1[i][j][k])+(Vt1[i][j+1][k]-Vt1[i][j][k]))/dex;
                        float PA3= PN[k+1]-P3;
                        float PA2= PN[k]-PN[k+1];
                        float PA1=PN[k]-P3;
                        Ws[i][j][0]=Ws[i][j][0]+VTS*(P4- PN[k]);
                        //**1
                        Wt[i][j][k]=(-V_TB*PA1+(Ws[i][j][0]*PA3/PA2-Wt[i][j][k+1]*PA2/PA3))/(PA3/PA2-PA2/PA3);
                        //**2
                        PHI[i][j][k]=(float) (PHI[i][j][k+1]+R/R0*TV*Math.pow((PN[k]/P0), RC)*(-PA2));
                    } else if (k==K-1){
                        float PA0= PN[k]-PN[k-1];
                        //***1
                        PHI[i][j][k]= (float) (PHI[i][j][k-1]+R/R0*TV*Math.pow(PN[k]/P0, RC)*PA0);
                        Wt[i][j][k]=0;
                    }else{
                        float PP1= PN[k-1]-PN[k];
                        float PP2= PN[k-1]-PN[k];
                        float PP3= PN[k]-PN[k+1];
                        PHI[i][j][k]=(float) ((R/R0*TV*Math.pow(PN[k]/P0, RC)*PP1+PHI[i][j][k-1]*PP3/PP2
                            -PHI[i][j][k+1]*PP2/PP3))/(PP3/PP2-PP2/PP3);
                        V_TB=-M*((Ut1[i+1][j][k]-Ut1[i][j][k])
                            +(Vt1[i][j+1][k]-Vt1[i][j][k])
                            +(Ut1[i+1][j][k+1]-Ut1[i][j][k+1])
                            +(Vt1[i][j+1][k+1]-Vt1[i][j][k+1]))/(2*dex);
                        Wt[i][j][k]=(-V_TB*PP1+Wt[i][j][k-1]*PP3/PP2-Wt[i][j][k+1]*PP2/PP1)
                            /(PP3/PP2- PP2/PP3);
                    }
                }
            }
        }
    }
    public static void printResult(){     
        for (int k = 0; k < K; k++) {
            for (int i = 0; i < X; i++) {
                for (int j = 0; j < Y; j++) {
                    System.out.println("U=");
                    System.out.println((i+1)+"\t"+(j+1)+"\t"+(k+1)+"\t"+Ut[i][j][k]);
                    System.out.println("V=");
                    System.out.println((i+1)+"\t"+(j+1)+"\t"+(k+1)+"\t"+Vt[i][j][k]);
                    System.out.println("TE=");
                    System.out.println((i+1)+"\t"+(j+1)+"\t"+(k+1)+"\t"+TE_t[i][j][k]);
                    System.out.println("Q=");
                    System.out.println((i+1)+"\t"+(j+1)+"\t"+(k+1)+"\t"+Q_t[i][j][k]);
                }
            }
        }       
    }
    
    public static void saveResult(){
        
    }
    public static void main(String[] args) {
        // TODO code application logic here
       
        System.out.println("CHUONG TRINH DU BAO VA CHUAN DOAN THOI TIET");
        System.out.println("-------------------------------------------");
        init();
        for(int n=0; n<N; n++){
            if (n==0) det1=det;
            //giai doan du doan
            DuDoan();
            System.out.println("Ket thuc viec du bao");
            //giai doan chuan doan
            ChuanDoan();
            System.out.println("Ket thuc viec chuan doan");
            System.out.println("Ket thuc tinh toan cho buoc thoi gian "+(n+1));
        }
        System.out.println("Ket thuc tinh toan");
        System.out.println("Lua chon cach in ket qua");
        System.out.println("1. Hien thi ra man hinh");
        System.out.println("2. Luu vao file ");
        System.out.println("Nhap lua chon: ");
        int result= scanner.nextInt();
        if (result==1) printResult();
        else if(result==2) saveResult();
    }
}
