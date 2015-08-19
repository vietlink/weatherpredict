
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package weatherpredict;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.math.BigDecimal;
import java.util.Scanner;

/**
 *
 * @author ngo
 */
public class WeatherPredict {

    /**
     * @param args the command line arguments
     */
    private static final int X=20;
    private static final int Y=20;
    private static final int K=5;
    private static final int N=50;
    private static final float R=287;
    private static final float CP=1005;
    private static final float M=1;
    
    private static final float[] PN= {1000, 900, 650, 400, 150};
    private static final int P0=1000;
    private static final float det=1200;
    private static final float dex=100000;
    private static final float F=(float) (14.584*Math.pow(10, -5));    
    private static final int NO_OF_DECIMAL=4;
    
    static float[][][] Ut= new float[K][X][Y];
    static float[][][] Ut1= new float[K][X][Y];
    static float[][][] Ut2= new float[K][X][Y];
    
    static float[][][] Vt= new float[K][X][Y];
    static float[][][] Vt1= new float[K][X][Y];
    static float[][][] Vt2= new float[K][X][Y];
    
    static float[][][] TE_t= new float[K][X][Y];
    static float[][][] TE_t1= new float[K][X][Y];
    static float[][][] TE_t2= new float[K][X][Y];
    static float[][][] Q_t= new float[K][X][Y];
    static float[][][] Q_t1= new float[K][X][Y];
    static float[][][] Q_t2= new float[K][X][Y];
    static float[][][] PHI= new float[K][X][Y];
    static float[][][] Wt= new float[K][X][Y];
    static float[][] PS= new float[X][Y];
    static float[][] PS_t1= new float[X][Y];
    static float[][] PS_t2= new float[X][Y];    
    static float[][] Ws= new float[X][Y];    
    static float fPi=(float) Math.PI;
    static float[] t= {18, 8, -41, -87, -233};
    static float[] Goc= {fPi/12, fPi/3, 2*fPi/3, 8*fPi/9, 35*fPi/36};
    static float[] Gio= {5f, 10.8f, 18f, 20.2f, 19.1f};
    static float[] Zcao= {110, 980, 4940, 7880, 14000};
    static float P1, P2, P3, P4, P5; 
    static float U1, U2, V1, V2;
    static float W1, W2, W9;
    static float R0,Ro, TV,V_TB;
    static float M1= M/(4*dex);
    static float RC= R/CP;
        
    static float det1;
    static Scanner scanner;   
    
    static String dataPath="C:\\Users\\ngo\\Desktop\\weatherpredict_result\\";
    public static void main(String[] args) throws InterruptedException {
        // TODO code application logic here
       
        System.out.println("CHUONG TRINH DU BAO VA CHUAN DOAN THOI TIET");
        System.out.println("-------------------------------------------");
        init();        
        for(int n=0; n<N; n++){            
            //giai doan du doan
            det1=2*det;
            if(n==0) {
                det1=det;
            }
            DuDoan();
//            System.out.println("Ket thuc viec du bao");            
            //giai doan chuan doan
            ChuanDoan();
//            System.out.println("Ket thuc viec chuan doan");
            System.out.println("Ket thuc tinh toan cho buoc thoi gian "+(n+1));
            saveResult(n);
       }
        System.out.println("Ket thuc tinh toan");        
//        printResult();        
    }
    
    public static void init() throws InterruptedException{
        float sgoc, cgoc, Tt, Ptb;
        System.out.println("Khoi tao so lieu");
        createFolder();
        for (int k = 0; k < K; k++) {
            sgoc=(float) Math.sin(Goc[k]);
//            roundValue(sgoc, NO_OF_DECIMAL);
            cgoc=(float) Math.cos(Goc[k]);
//            roundValue(cgoc, NO_OF_DECIMAL);
            Tt= 273+t[k];
            if(k==(K-1)) Ptb=25;
            else Ptb=(PN[k]+PN[k+1])/2;
            for (int i = 0; i < X; i++) {
                for (int j = 0; j < Y; j++) {
                    Ut[k][i][j]= Gio[k]*cgoc;
//                    Ut[k][i][j]=roundValue(Ut[k][i][j], NO_OF_DECIMAL);
                    Vt[k][i][j]= Gio[k]*sgoc;
//                    Vt[k][i][j]=roundValue(Vt[k][i][j], NO_OF_DECIMAL);
                    TE_t[k][i][j]= (float) (Tt*Math.pow(1000/Ptb, 0.286));
//                    TE_t[k][i][j]=roundValue(TE_t[k][i][j], NO_OF_DECIMAL);
                    Q_t[k][i][j]= 0;
                    PHI[k][i][j]=Zcao[k];
                    Wt[k][i][j]=0;
                    Ut1[k][i][j]=Ut[k][i][j]; 
                    Vt1[k][i][j]=Vt[k][i][j];
                    TE_t1[k][i][j]=TE_t[k][i][j]; 
                    Q_t1[k][i][j]=Q_t[k][i][j];
                }
            }
        }
        for (int i = 0; i < X; i++) {
            for (int j = 0; j < Y; j++) {
                PS[i][j]=P0;
                Ws[i][j]=0;
                PS_t1[i][j]=PS[i][j];
            }
        }
        saveFile3D("init\\U.txt", Ut);
        saveFile3D("init\\V.txt", Vt);
        saveFile3D("init\\TE.txt", TE_t);
        saveFile3D("init\\Q.txt", Q_t);
        saveFile3D("init\\PHI.txt", PHI);
        saveFile3D("init\\Wt.txt", Wt);
        saveFile2D("init\\Ws.txt", Ws);
        saveFile2D("init\\PS.txt", PS);
//        print3D(K, X, Y, Ut1);
//        print3D(K, X, Y, Vt1);
//        print3D(K, X, Y, TE_t1);
//        print3D(K, X, Y, Q_t1);
//        print3D(K, X, Y, PHI);
//        print3D(K, X, Y, Wt);
//        print2D(X, Y, Ws);
//        print2D(X, Y, PS);
        System.out.println("Khoi tao thanh cong");
        
    }
    public static void createFolder(){
        for (int i = 0; i < N; i++) {
            new File(dataPath+"tinhtoan\\"+(i+1)).mkdir();            
            }
    }
    public static void DuDoan(){       
        for(int k=0; k<K; k++){                
            for (int i=1;i<X-1; i++){
                for(int j=1; j<Y-1; j++){
                    //2.1                    
                    Ut2[k][i][j]=-M*(PHI[k][i+1][j]-PHI[k][i][j])/dex
                        -M1*((Ut1[k][i][j]+Ut1[k][i+1][j])*(Ut1[k][i+1][j]-Ut1[k][i][j])+(Ut1[k][i][j]+Ut1[k][i-1][j])*(Ut1[k][i][j]-Ut1[k][i-1][j]))
                        +F*(Vt1[k][i][j]+Vt1[k][i][j+1]+Vt1[k][i-1][j]+Vt1[k][i-1][j+1])/4
                        -M1*((Vt1[k][i][j]+Vt1[k][i+1][j])*(Ut1[k][i][j+1]-Ut1[k][i][j])+(Vt1[k][i][j-1]+Vt1[k][i+1][j-1])*(Ut1[k][i][j-1]-Ut1[k][i][j]));  
//                    Ut2[k][i][j]=roundValue(Ut2[k][i][j], NO_OF_DECIMAL);
                    //2.2
                    Vt2[k][i][j]= -M*(PHI[k][i][j+1]-PHI[k][i][j])/dex
                        - M1*((Ut1[k][i][j]+Ut1[k][i][j+1])*(Vt1[k][i+1][j]-Vt1[k][i][j])+(Ut1[k][i-1][j]+Ut1[k][i-1][j+1])*(Vt1[k][i][j]-Vt1[k][i-1][j]))
                        - M1*((Vt1[k][i][j]+Vt1[k][i][j+1])*(Vt1[k][i][j+1]-Vt1[k][i][j])+(Vt1[k][i][j]+Vt1[k][i][j-1])*(Vt1[k][i][j]-Vt1[k][i][j-1]))
                        -F*(Ut1[k][i][j]+Ut1[k][i][j+1]+Ut1[k][i][j-1]+Vt1[k][i+1][j-1])/4;      
//                    Vt2[k][i][j]=roundValue(Vt2[k][i][j], NO_OF_DECIMAL);
                    if(k==K-1){
                        //***
                        P1= (PN[k]+PN[k-1])/2;
//                        P1=roundValue(P1, NO_OF_DECIMAL);
                        P2=25;  
                        P5=(PN[k]-P2)/(2*(P1-P2));
//                        P5=roundValue(P5, NO_OF_DECIMAL);
                        //***1
                        Ut2[k][i][j]= Ut2[k][i][j]
                                -((Wt[k-1][i+1][j]+Wt[k-1][i][j])*P5*(Ut1[k][i][j]-Ut1[k-1][i][j]))/(PN[k]-PN[k-1]);
//                        Ut2[k][i][j]= roundValue(Ut2[k][i][j], NO_OF_DECIMAL);
                        //***2
                        Vt2[k][i][j]= Vt2[k][i][j]
                                -((Wt[k-1][i][j+1]+Wt[k-1][i][j])*P5*(Vt1[k][i][j]-Vt1[k-1][i][j]))/(PN[k]-PN[k-1]);
//                        Vt2[k][i][j]= roundValue(Vt2[k][i][j], NO_OF_DECIMAL);
                        //***3
                        TE_t2[k][i][j]=0;
                        Q_t2[k][i][j]=0;
                    } else{
                        U1= (Ut1[k][i+1][j]+Ut1[k+1][i+1][j])/2;
//                        U1=roundValue(U1, NO_OF_DECIMAL);
                        U2= (Ut1[k][i][j]+Ut1[k+1][i][j])/2;
//                        U2=roundValue(U2, NO_OF_DECIMAL);
                        V1= (Vt1[k][i][j+1]+Vt1[k+1][i][j+1])/2;
//                        V1=roundValue(V1, NO_OF_DECIMAL);
                        V2=(Vt1[k][i][j]+Vt1[k+1][i][j])/2;
//                        V2=roundValue(V2, NO_OF_DECIMAL);
                        //2.3
                        TE_t2[k][i][j]=-2*M1*(U1*(TE_t1[k][i+1][j]-TE_t1[k][i][j])
                                -U2*(TE_t1[k][i][j]-TE_t1[k][i-1][j])
                                +V1*(TE_t1[k][i][j+1]-TE_t1[k][i][j])
                                -V2*(TE_t1[k][i][j]-TE_t1[k][i][j-1]));
//                        TE_t2[k][i][j]=roundValue(TE_t2[k][i][j], NO_OF_DECIMAL);
                        //2.4
                        Q_t2[k][i][j]=-2*M1*(U1*(Q_t1[k][i+1][j]-Q_t1[k][i][j])
                                -U2*(Q_t1[k][i][j]-Q_t1[k][i-1][j])
                                +V1*(Q_t1[k][i][j+1]-Q_t1[k][i][j])
                                -V2*(Q_t1[k][i][j]-Q_t1[k][i][j-1]));
//                        Q_t2[k][i][j]=roundValue(Q_t2[k][i][j], NO_OF_DECIMAL);
                        if(k==0){
                        P1= (PN[k+1]+PN[k])/2;
//                        P1=roundValue(P1, NO_OF_DECIMAL);
                        P2= (PN[k+2]+PN[k+1])/2;    
//                        P2=roundValue(P2, NO_OF_DECIMAL);
                        W9= (Wt[k+1][i][j]*P1+Wt[k][i][j]*P2)/(P1+P2);
//                        W9=roundValue(W9, NO_OF_DECIMAL);
                        W1= ((Wt[k][i+1][j]+Wt[k][i][j])+(Ws[i][j]+Ws[i+1][j]))/4;
//                        W1=roundValue(W1, NO_OF_DECIMAL);
                        W2=((Wt[k][i][j+1]+Wt[k][i][j])+(Ws[i][j+1]+Ws[i][j]))/4;
//                        W2=roundValue(W2, NO_OF_DECIMAL);
                        //**1
                        Ut2[k][i][j]=Ut2[k][i][j]-W1*(Ut1[k+1][i][j]-Ut1[k][i][j])/(PN[k+1]-PN[k]);
//                        Ut2[k][i][j]=roundValue(Ut2[k][i][j], NO_OF_DECIMAL);
                        System.err.println("Ut2["+k+"]["+i+"]["+j+"]= "+Ut2[k][i][j]);
                        //**2
                        Vt2[k][i][j]=Vt2[k][i][j]-W2*(Vt1[k+1][i][j]-Vt1[k][i][j])/(PN[k+1]-PN[k]);
//                        System.err.println("Vt2["+k+"]["+i+"]["+j+"]= "+Vt2[k][i][j]);
//                        Vt2[k][i][j]=roundValue(Vt2[k][i][j], NO_OF_DECIMAL);
                        //**3
                        TE_t2[k][i][j]=TE_t2[k][i][j]-Wt[k][i][j]*(TE_t1[k+1][i][j]-TE_t1[k][i][j])/(PN[k+1]-PN[k]);
//                        System.err.println("TE_t2["+k+"]["+i+"]["+j+"]= "+TE_t2[k][i][j]);
//                        TE_t2[k][i][j]=roundValue(TE_t2[k][i][j], NO_OF_DECIMAL);
                        //**4
                        Q_t2[k][i][j]= Q_t2[k][i][j]-Wt[k][i][j]*(Q_t1[k+1][i][j]-Q_t1[k][i][j])/(PN[k+1]-PN[k]);
//                        System.err.println("Q_t2["+k+"]["+i+"]["+j+"]= "+Q_t2[k][i][j]);
//                        Q_t2[k][i][j]= roundValue(Q_t2[k][i][j], NO_OF_DECIMAL);
                        //**5
                        PS_t2[i][j]=-2*M1*(Ut1[k][i+1][j]*(PS_t1[i+1][j]-PS_t1[i][j])
                                -Ut1[k][i][j]*(PS_t1[i][j]-PS_t1[i-1][j])
                                +Vt1[k][i][j+1]*(PS_t1[i][j+1]-PS_t1[i][j])
                                -Vt1[k][i][j]*(PS_t1[i][j]-PS_t1[i][j-1]));
//                        PS_t2[i][j]=roundValue(PS_t2[i][j], NO_OF_DECIMAL);
                        PS_t2[i][j]=PS_t2[i][j]+W9-(((Ut1[k+1][i+1][j]-Ut1[k+1][i][j])
                                +(Vt1[k+1][i][j+1]-Vt1[k+1][i][j])+(Ut1[k][i+1][j]-Ut1[k][i][j])
                                +(Vt1[k][i][j+1]-Vt1[k][i][j]))*2*M1)*(PS_t1[i][j]-900);
//                        PS_t2[i][j]=roundValue(PS_t2[i][j], NO_OF_DECIMAL);
//                       System.err.println("PS_t2["+i+"]["+j+"]= "+PS_t2[i][j]);
                        } else{                    
                        P1= (PN[k-1]+PN[k])/2;
                        P2= (PN[k]+PN[k+1])/2;
                        //2.5
                        Ut2[k][i][j]= Ut2[k][i][j]-((Wt[k][i+1][j]+Wt[k][i][j])*(PN[k]-P2)
                                +(Wt[k+1][i+1][j]+Wt[k+1][i][j])*(P1-PN[k])/(2*(P1-P2))
                                *(Ut1[k-1][i][j]-Ut1[k][i][j])*(PN[k]-PN[k+1])/(PN[k-1]-PN[k])
                                +(Ut1[k][i][j]-Ut1[k+1][i][j])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1]))/(PN[k-1]-PN[k+1]);
//                        Ut2[k][i][j]= roundValue(Ut2[k][i][j], NO_OF_DECIMAL);
                        //2.6
                        Vt2[k][i][j]= Vt2[k][i][j]-((Wt[k][i][j+1]+Wt[k][i][j])*(PN[k]-P2)
                                +(Wt[k+1][i][j+1]+Wt[k+1][i][j])*(P1-PN[k])/(2*(P1-P2))
                                *(Vt1[k-1][i][j]-Vt1[k][i][j])*(PN[k]-PN[k+1])/(PN[k-1]-PN[k])
                                +(Vt1[k][i][j]-Vt1[k+1][i][j])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1]))/(PN[k-1]-PN[k+1]);
//                        Vt2[k][i][j]= roundValue(Vt2[k][i][j], NO_OF_DECIMAL);
                        //2.7
                        TE_t2[k][i][j]=TE_t2[k][i][j]-Wt[k][i][j]*((TE_t1[k-1][i][j]-TE_t1[k][i][j])*(PN[k]-PN[k+1])/(PN[k-1]-PN[k])
                                +(TE_t1[k][i][j]-TE_t[k+1][i][j])*(PN[k+1]-PN[k])/(PN[k]-PN[k+1]))/(PN[k-1]-PN[k+1]);
//                        TE_t2[k][i][j]=roundValue(TE_t2[k][i][j], NO_OF_DECIMAL);
                        //2.8
                        Q_t2[k][i][j]= Q_t2[k][i][j]-Wt[k][i][j]*((Q_t1[k-1][i][j]-Q_t1[k][i][j])*(PN[k]-PN[k+1])/(PN[k-1]-PN[k])
                                +(Q_t1[k][i][j]-Q_t1[k+1][i][j])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1]))/(PN[k-1]-PN[k+1]);
//                        Q_t2[k][i][j]=roundValue(Q_t2[k][i][j], NO_OF_DECIMAL);
                        }
                    }                                                                                          
                }
            }
        }       
        //xu li ket qua du doan
        for (int k = 0; k < K; k++) {
            for (int i = 1; i < X-1; i++) {
                for (int j = 1; j < Y-1; j++) {
                    Ut2[k][i][j]= Ut[k][i][j]+det1*Ut2[k][i][j];
//                    Ut2[k][i][j]=roundValue(Ut2[k][i][j], NO_OF_DECIMAL);
                    Vt2[k][i][j]= Vt[k][i][j]+det1*Vt2[k][i][j];
//                    Vt2[k][i][j]= roundValue(Vt2[k][i][j], NO_OF_DECIMAL);
                    TE_t2[k][i][j]= TE_t[k][i][j]+det1*TE_t2[k][i][j];
//                    TE_t2[k][i][j]=roundValue(TE_t2[k][i][j], NO_OF_DECIMAL);
                    Q_t2[k][i][j]= Q_t[k][i][j]+det1*Q_t2[k][i][j];                        
//                    Q_t2[k][i][j]=roundValue(Q_t2[k][i][j], NO_OF_DECIMAL);
                    Ut[k][i][j]=Ut1[k][i][j]; Ut1[k][i][j]=Ut2[k][i][j];
                    Vt[k][i][j]=Vt1[k][i][j]; Vt1[k][i][j]=Vt2[k][i][j];
                    TE_t[k][i][j]=TE_t1[k][i][j]; TE_t1[k][i][j]=TE_t2[k][i][j];
                    Q_t[k][i][j]=Q_t1[k][i][j]; Q_t1[k][i][j]=Q_t2[k][i][j];
                }
            }
        }
        for (int i = 1; i < X-1; i++) {
            for (int j = 1; j < Y-1; j++) {
                PS_t2[i][j]= PS[i][j]+det1*PS_t2[i][j];
//                PS_t2[i][j]=roundValue(PS_t2[i][j], NO_OF_DECIMAL);
                PS[i][j]=PS_t1[i][j];
                PS_t1[i][j]=PS_t2[i][j];
            }                
        }
    }
    
    public static void ChuanDoan(){
        for(int k=0; k<K; k++){
            for (int i=1; i<X-1; i++){
                for(int j=1; j<Y-1; j++){
                    R0= (float) (TE_t1[k][i][j]*Math.pow(PN[k]/P0, RC));
                    Ro= PN[k]/(R*R0);
                    TV=(float) (TE_t1[k][i][j]*(1+0.61*Q_t1[k][i][j]));
                    if(k==0){
                        P1=(PN[k]-PN[k+1])/2;
                        P2=(PN[k+1]-PN[k+2])/2;
                        P3=(PN[k+1]+PN[k+2])/2;                        
                            //**
                        V_TB=-M*(((Ut1[k][i+1][j]-Ut1[k][i][j])
                            +(Vt1[k][i][j+1]-Vt1[k][i][j]))*P2     
                            +((Ut1[k+1][i+1][j]-Ut1[k+1][i][j])
                            +(Vt1[k+1][i][j+1]-Vt1[k+1][i][j]))*P1)/(dex*(P1+P2));          
//                        V_TB=roundValue(V_TB, NO_OF_DECIMAL);
                        float PA3= PN[k+1]-P3;
                        float PA2= PN[k]-PN[k+1];
                        float PA1=PN[k]-P3;                        
                        //**1
                        Wt[k][i][j]=(-V_TB*PA1+(Ws[i][j]*PA3/PA2-Wt[k+1][i][j]*PA2/PA3))/(PA3/PA2-PA2/PA3);
//                        Wt[k][i][j]=roundValue(Wt[k][i][j], NO_OF_DECIMAL);
                        //**2
                        PHI[k][i][j]=(float) (PHI[k+1][i][j]+R/Ro*TV*Math.pow((PN[k]/P0), RC)*(-PA2));
//                        PHI[k][i][j]=roundValue(PHI[k][i][j], NO_OF_DECIMAL);
                    } else if (k==K-1){
                        float PA0= PN[k]-PN[k-1];
                        //***1
                        PHI[k][i][j]=(float) (PHI[k-1][i][j]+R/Ro*TV*Math.pow(PN[k]/P0, RC)*PA0);
//                        PHI[k][i][j]=roundValue(PHI[k][i][j], NO_OF_DECIMAL);
                        Wt[k][i][j]=0;
                    }else{
                        float PP1= PN[k-1]-PN[k+1];
                        float PP2= PN[k-1]-PN[k];
                        float PP3= PN[k]-PN[k+1];
                        PHI[k][i][j]=(float) (((R/Ro*TV*Math.pow(PN[k]/P0, RC)*PP1+PHI[k-1][i][j]*PP3/PP2
                                -PHI[k+1][i][j]*PP2/PP3))/(PP3/PP2-PP2/PP3));
//                        PHI[k][i][j]=roundValue(PHI[k][i][j], NO_OF_DECIMAL);
                        V_TB=-M*((Ut1[k][i+1][j]-Ut1[k][i][j])
                            +(Vt1[k][i][j+1]-Vt1[k][i][j])
                            +(Ut1[k+1][i+1][j]-Ut1[k+1][i][j])
                            +(Vt1[k+1][i][j+1]-Vt1[k+1][i][j]))/(2*dex);
//                        V_TB=roundValue(V_TB, NO_OF_DECIMAL);
                        Wt[k][i][j]=(-V_TB*PP1+Wt[k-1][i][j]*PP3/PP2-Wt[k+1][i][j]*PP2/PP1)
                            /(PP3/PP2- PP2/PP3);
//                        Wt[k][i][j]=roundValue(Wt[k][i][j], NO_OF_DECIMAL);
                    }
                }
            }
        }
        P4=(PN[1]+PN[2])/2;
        for (int i = 1; i < X-1; i++) {
            for (int j = 1; j < Y-1; j++) {
                float VTS= M*((Ut1[1][i+1][j]-Ut1[1][i][j])+(Vt1[1][i][j+1]-Vt1[1][i][j]))/dex;
//                VTS=roundValue(VTS, NO_OF_DECIMAL);
                Ws[i][j]=Wt[1][i][j]+VTS*(P4- PN[1]);
//                Ws[i][j]=roundValue(Ws[i][j], NO_OF_DECIMAL);
            }
        }
    }
    public static void printResult(){   
        System.err.println("U=");
        print3D(K, X, Y, Ut2);
        
        System.err.println("V=");
        
        print3D(K, X, Y, Vt2);
        System.err.println("TE=");
        
        print3D(K, X, Y, TE_t2);
        System.err.println("Q=");
        
        print3D(K, X, Y, Q_t2);
        System.err.println("Wt=");
        
        print3D(K, X, Y, Wt);
        System.err.println("PHI=");
        
        print3D(K, X, Y, PHI);
        
        System.err.println("PS=");
        print2D(X, Y, PS_t2);
        
        System.err.println("Ws=");
        
        print2D(X, Y, Ws);        
    }
    public static void saveFile3D(String fileName,float[][][] array ){
        try {            
            BufferedWriter bw= null;
            File file= new File(dataPath+fileName);
            if(!file.exists()){
                file.createNewFile();
            }
            bw= new BufferedWriter(new FileWriter(file));
            for (int k = 0; k < array.length; k++) {
                bw.write("k= "+(k+1));
                bw.newLine();
                for (int i = 0; i < array[1].length; i++) {
                    for (int j = 0; j < array[0][1].length; j++) {
                        bw.write((i+1)+"\t"+(j+1)+"\t"+array[k][i][j]);
                        bw.newLine();
//                        bw.write(array[k][i][j]+" ");                        
                    }
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
        }
    }
     public static void saveFile2D(String fileName,float[][] array ){
        try {            
            BufferedWriter bw= null;
            File file= new File(dataPath+fileName);
            if(!file.exists()){
                file.createNewFile();
            }
            bw= new BufferedWriter(new FileWriter(file));
            for (int i = 0; i < array.length; i++) {
                for (int j = 0; j < array[1].length; j++) {                    
//                        bw.write(array[i][j]+" ");                                            
                    bw.write((i+1)+"\t"+(j+1)+"\t"+array[i][j]);
                    bw.newLine();
                }                
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
        }
    }
    public static void print3D(int k, int i, int j, float[][][] array){
        for (int l = 0; l < k; l++) {
            System.out.println("k= "+(l+1));
            print2D(i, j, array[l]);
        }
    }
    public static void print2D(int i, int j, float[][] array){
        for (int k = 0; k < i; k++) {
            for (int l = 0; l < j; l++) {
//                System.out.println((k+1)+"\t"+(l+1)+"\t"+array[k][l]);
                System.out.print(BigDecimal.valueOf(array[k][l])+" ");
            }
            System.out.println();
        }
    }
    //luu gia tri tai tung buoc n vao cac folder
    public static void saveResult(int n){
        saveFile3D("tinhtoan\\"+n+"\\Ut2_"+N+".txt", Ut2);
        saveFile3D("tinhtoan\\"+n+"\\Vt2_"+N+".txt", Vt2);
        saveFile3D("tinhtoan\\"+n+"\\TE_t2_"+N+".txt", TE_t2);
        saveFile3D("tinhtoan\\"+n+"\\Q_t2_"+N+".txt", Q_t2);
        saveFile3D("tinhtoan\\"+n+"\\Wt_"+N+".txt", Wt);
        saveFile3D("tinhtoan\\"+n+"\\PHI_"+N+".txt", PHI);
        saveFile2D("tinhtoan\\"+n+"\\PS_t2_"+N+".txt", PS_t2);
        saveFile2D("tinhtoan\\"+n+"\\Ws_"+N+".txt", Ws);
    }
    public static float roundValue(float value, int n){
        float decimal_place= (float) Math.pow(10, n);
        return (Math.round(value*decimal_place))/decimal_place;
    }
    
}
