
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package weatherpredict;

import com.sun.org.apache.xml.internal.serialize.LineSeparator;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
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
    private static final double R=287;
    private static final double CP=1005;
    private static final double M=1;
    private static final double N=50;
    private static final double[] PN= {1000, 900, 650, 400, 150};
    private static final int P0=1000;
    private static final double det=1200;
    private static final double dex=100000;
    private static final double F=(double) (14.584*Math.pow(10, -5));    
    
    static double[][][] Ut= new double[K][X][Y];
    static double[][][] Ut1= new double[K][X][Y];
    static double[][][] Ut2= new double[K][X][Y];
    
    static double[][][] Vt= new double[K][X][Y];
    static double[][][] Vt1= new double[K][X][Y];
    static double[][][] Vt2= new double[K][X][Y];
    
    static double[][][] TE_t= new double[K][X][Y];
    static double[][][] TE_t1= new double[K][X][Y];
    static double[][][] TE_t2= new double[K][X][Y];
    static double[][][] Q_t= new double[K][X][Y];
    static double[][][] Q_t1= new double[K][X][Y];
    static double[][][] Q_t2= new double[K][X][Y];
    static double[][][] PHI= new double[K][X][Y];
    static double[][][] Wt= new double[K][X][Y];
    static double[][] PS= new double[X][Y];
    static double[][] PS_t1= new double[X][Y];
    static double[][] PS_t2= new double[X][Y];    
    static double[][] Ws= new double[X][Y];    
    
    static double[] t= {18, 8, -41, -87, -233};
    static double[] Goc= {Math.PI/12, Math.PI/3, 2*Math.PI/3, 8*Math.PI/9, 35*Math.PI/36};
    static double[] Gio= {5, 10.8, 18, 20.2, 19.1};
    static double[] Zcao= {110, 980, 4940, 7880, 14000};
    static double P1, P2, P3, P4, P5; 
    static double U1, U2, V1, V2;
    static double W1, W2, W9;
    static double R0,Ro, TV,V_TB;
    static double M1= M/(4*dex);
    static double RC= R/CP;
        
    static double det1;
    static Scanner scanner;   
    
    static String dataPath="C:\\Users\\ngo\\Desktop\\temp\\";
    public static void init(){
        double sgoc, cgoc, Tt, Ptb;
        System.out.println("Khoi tao so lieu");
        for (int k = 0; k < K; k++) {
            sgoc=Math.sin(Goc[k]);
            roundValue(sgoc, 3);
            cgoc=Math.cos(Goc[k]);
            roundValue(cgoc, 3);
            Tt= 273+t[k];
            if(k==(K-1)) Ptb=25;
            else Ptb=(PN[k]+PN[k+1])/2;
            for (int i = 0; i < X; i++) {
                for (int j = 0; j < Y; j++) {
                    Ut[k][i][j]= Gio[k]*cgoc;
                    Ut[k][i][j]=roundValue(Ut[k][i][j], 3);
                    Vt[k][i][j]= Gio[k]*sgoc;
                    Vt[k][i][j]=roundValue(Vt[k][i][j], 3);
                    TE_t[k][i][j]= Tt*Math.pow(1000/Ptb, 0.286);
                    TE_t[k][i][j]=roundValue(TE_t[k][i][j], 3);
                    Q_t[k][i][j]= 0;
                    PHI[k][i][j]=Zcao[k];
                    Wt[k][i][j]=0;
                    Ut1[k][i][j]=Ut[k][i][j]; Vt1[k][i][j]=Vt[k][i][j];
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
        System.out.println("Khoi tao thanh cong");
        System.out.println("Ut");
//        print3D(K, X, Y, Ut1);
//        print3D(K, X, Y, Vt);
//        print3D(K, X, Y, TE_t);
        print3D(K, X, Y, Q_t);
    }
      
    public static void DuDoan(){
        System.out.println(det1);
        for(int k=0; k<K; k++){                
            for (int i=1;i<X-1; i++){
                for(int j=1; j<Y-1; j++){
                    //2.1                    
                    Ut2[k][i][j]=-M*(PHI[k][i+1][j]-PHI[k][i][j])/dex
                        -M1*((Ut1[k][i][j]+Ut1[k][i+1][j])*(Ut1[k][i+1][j]-Ut1[k][i][j])+(Ut1[k][i][j]+Ut1[k][i-1][j])*(Ut1[k][i][j]-Ut1[k][i-1][j]))
                        +F*(Vt1[k][i][j]+Vt1[k][i][j-1]+Vt1[k][i+1][j]+Vt1[k][i+1][j-1])/4
                        -M1*((Vt1[k][i][j]+Vt1[k][i+1][j])*(Ut1[k][i][j+1]-Ut1[k][i][j])+(Vt1[k][i][j-1]+Vt1[k][i+1][j-1])*(Ut1[k][i][j-1]-Ut1[k][i][j]));  
                    Ut2[k][i][j]=roundValue(Ut2[k][i][j], 3);
                    //2.2
                    Vt2[k][i][j]= -M*(PHI[k][i][j+1]-PHI[k][i][j])/dex
                        - M1*((Ut1[k][i][j]+Ut1[k][i][j+1])*(Vt1[k][i+1][j]-Vt1[k][i][j])+(Ut1[k][i-1][j]+Ut1[k][i-1][j+1])*(Vt1[k][i][j]-Vt1[k][i-1][j]))
                        - M1*((Vt1[k][i][j]+Vt1[k][i][j+1])*(Vt1[k][i][j+1]-Vt1[k][i][j])+(Vt1[k][i][j]+Vt1[k][i][j-1])*(Vt1[k][i][j]-Vt1[k][i][j-1]))
                        -F*(Ut1[k][i][j]+Ut1[k][i-1][j]+Ut1[k][i][j+1]+Vt1[k][i-1][j+1])/4;      
                    Vt2[k][i][j]=roundValue(Vt2[k][i][j], 3);                    
                     if(k==K-1){
                        //***
                        P1= (PN[k]+PN[k-1])/2;                            
                        P2=25;  
                        P5=(PN[k]-P2)/(2*(P1-P2));
                        //***1
                        Ut2[k][i][j]= Ut2[k][i][j]
                                -((Wt[k-1][i+1][j]+Wt[k-1][i][j])*P5*(Ut1[k][i][j]-Ut1[k-1][i][j]))/(PN[k]-PN[k-1]);
                        Ut2[k][i][j]= roundValue(Ut2[k][i][j], 3);
                        //***2
                        Vt2[k][i][j]= Vt2[k][i][j]
                                -((Wt[k-1][i][j+1]+Wt[k-1][i][j])*P5*(Vt1[k][i][j]-Vt1[k-1][i][j]))/(PN[k]-PN[k-1]);
                        Vt2[k][i][j]= roundValue(Vt2[k][i][j], 3);
                        //***3
                        TE_t2[k][i][j]=0;
                        Q_t2[k][i][j]=0;
                    }                    
                     else if(k==0){
                        P1= (PN[k+1]+PN[k])/2;                            
                        P2= (PN[k+2]+PN[k+1])/2;    
                        W9= (Wt[k+2][i][j]*P1+Wt[k][i][j]*P2)/(P1+P2);
                        W1= ((Wt[k][i+1][j]+Wt[k][i][j])+(Ws[i][j]+Ws[i+1][j]))/4;
                        W2=((Wt[k][i][j+1]+Wt[k][i][j])+(Ws[i][j+1]+Ws[i][j]))/4;
                        //**1
                        Ut2[k][i][j]=Ut2[k][i][j]-W1*(Ut1[k+1][i][j]-Ut1[k][i][j])/(PN[k+1]-PN[k]);
                        Ut2[k][i][j]=roundValue(Ut2[k][i][j], 3);
                        //**2
                        Vt2[k][i][j]=Vt2[k][i][j]-W2*(Vt1[k+1][i][j]-Vt1[k][i][j])/(PN[k+1]-PN[k]);
                        Vt2[k][i][j]=roundValue(Vt2[k][i][j], 3);
                        //**3
                        TE_t2[k][i][j]=TE_t2[k][i][j]-Wt[k][i][j]*(TE_t1[k+1][i][j]-TE_t1[k][i][j])/(PN[k+1]-PN[k]);
                        TE_t2[k][i][j]=roundValue(TE_t2[k][i][j], 3);
                        //**4
                        Q_t2[k][i][j]= Q_t2[k][i][j]-Wt[k][i][j]*(Q_t1[k+1][i][j]-Q_t1[k][i][j])/(PN[k+1]-PN[k]);
                        Q_t2[k][i][j]= roundValue(Q_t2[k][i][j], 3);
                        //**5
                        PS_t2[i][j]=-2*M1*(Ut1[k][i+1][j]*(PS_t1[i+1][j]-PS_t1[i][j])
                                -Ut1[k][i][j]*(PS_t1[i][j]-PS_t1[i-1][j])
                                -Vt1[k][i][j]*(PS_t1[i][j]-PS_t1[i][j-1]));
                        PS_t2[i][j]=roundValue(PS_t2[i][j], 3);
                        PS_t2[i][j]=PS_t2[i][j]+W9-(M*((Ut1[k+1][i+1][j]-Ut1[k][i][j])
                                +(Vt1[k+1][i][j+1]-Vt1[k+1][i][j])+(Ut1[k][i+1][j]-Ut1[k][i][j])
                                +(Vt1[k][i][j+1]-Vt1[k][i][j]))*2*M1)*(PS_t1[i][j]-900);
                        PS_t2[i][j]=roundValue(PS_t2[i][j], 3);
                    }                    
                    else{
                        U1= (Ut1[k][i+1][j]+Ut1[k+1][i+1][j])/2;
                        U2= (Ut1[k][i][j]+Ut1[k+1][i][j])/2;
                        V1= (Vt1[k][i][j+1]+Vt1[k+1][i][j+1])/2;
                        V2=(Vt1[k][i][j]+Vt1[k+1][i][j])/2;
                        //2.3
                        TE_t2[k][i][j]=-2*M1*(U1*(TE_t1[k][i+1][j]-TE_t1[k][i][j])
                                -U2*(TE_t1[k][i][j]-TE_t1[k][i-1][j])
                                +V1*(TE_t1[k][i][j+1]-TE_t1[k][i][j])
                                -V2*(TE_t1[k][i][j]-TE_t1[k][i][j-1]));
                        TE_t2[k][i][j]=roundValue(TE_t2[k][i][j], 3);
    //                    System.out.println("2.3");
                        //2.4
                        Q_t2[k][i][j]=-2*M1*(U1*(Q_t1[k][i+1][j]-Q_t1[k][i][j])
                                -U2*(Q_t1[k][i][j]-Q_t1[k][i-1][j])
                                +V1*(Q_t1[k][i][j+1]-Q_t1[k][i][j])
                                -V2*(Q_t1[k][i][j]-Q_t1[k][i][j-1]));
                        Q_t2[k][i][j]=roundValue(Q_t2[k][i][j], 3);
//                    System.out.println("2.4");
                        P1= (PN[k-1]+PN[k])/2;
                        P2= (PN[k]+PN[k+1])/2;
                        //2.5
                        Ut2[k][i][j]= Ut2[k][i][j]-(((Wt[k][i+1][j]+Wt[k][i][j])*(PN[k]-P2)
                                +(Wt[k+1][i+1][j]+Wt[k+1][i][j])*(P1-PN[k]))/(2*(P1-P2))
                                *((Ut1[k-1][i][j]-Ut1[k][i][j])*(PN[k]-PN[k+1])/(PN[k-1]-PN[k]))
                                +((Ut1[k][i][j]-Ut1[k+1][i][j])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1])))/(PN[k-1]-PN[k+1]);
                        Ut2[k][i][j]= roundValue(Ut2[k][i][j], 3);
                        //2.6
                        Vt2[k][i][j]= Vt2[k][i][j]-(((Wt[k][i][j+1]+Wt[k][i][j])*(PN[k]-P2)
                                +(Wt[k+1][i][j+1]+Wt[k+1][i][j])*(P1-PN[k]))/(2*(P1-P2))
                                *((Vt1[k-1][i][j]-Vt1[k+1][i][j])*(PN[k]-PN[k+1])/(PN[k-1]-PN[k]))
                                +((Vt1[k][i][j]-Vt1[k+1][i][j])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1])))/(PN[k-1]-PN[k+1]);
                        Vt2[k][i][j]= roundValue(Vt2[k][i][j], 3);
                        //2.7
                        TE_t2[k][i][j]=TE_t2[k][i][j]-Wt[k][i][j]*((TE_t1[k-1][i][j]-TE_t1[k][i][j])
                                *(PN[k]-PN[k+1])/(PN[k-1]-PN[k])+(TE_t1[k][i][j]-TE_t[k+1][i][j])*(PN[k-1]-PN[k])/(PN[k]-PN[k+1]))/(PN[k-1]-PN[k+1]);
                        TE_t2[k][i][j]=roundValue(TE_t2[k][i][j], 3);
                        //2.8
                        Q_t2[k][i][j]= Q_t2[k][i][j]-Wt[k][i][j]*((Q_t1[k-1][i][j]-Q_t1[k][i][j])*
                                (PN[k]-PN[k+1])/(PN[k-1]-PN[k])+(Q_t1[k][i][j]-Q_t1[k+1][i][j])*
                                (PN[k-1]-PN[k])/(PN[k]-PN[k+1]))/(PN[k-1]-PN[k+1]);
                        Q_t2[k][i][j]=roundValue(Q_t2[k][i][j], 3);
                    }
                    
                }
            }
        }       
        //xu li ket qua du doan
        for (int k = 0; k < K; k++) {
            for (int i = 1; i < X-1; i++) {
                for (int j = 1; j < Y-1; j++) {
                    Ut2[k][i][j]= Ut[k][i][j]+det1*Ut2[k][i][j];
                    Vt2[k][i][j]= Vt[k][i][j]+det1*Vt2[k][i][j];
                    TE_t2[k][i][j]= TE_t[k][i][j]+det1*TE_t2[k][i][j];
                    Q_t2[k][i][j]= Q_t[k][i][j]+det1*Q_t2[k][i][j];                        
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
                PS[i][j]=PS_t1[i][j];
                PS_t1[i][j]=PS_t2[i][j];
            }                
        }
    }
    
    public static void ChuanDoan(){
        for(int k=0; k<K; k++){
            for (int i=1; i<X-1; i++){
                for(int j=1; j<Y-1; j++){
                    R0= (double) (TE_t1[k][i][j]*Math.pow(PN[k]/P0, RC));
                    Ro= PN[k]/(R*R0);
                    TV=(double) (TE_t1[k][i][j]*(1+0.61*Q_t1[k][i][j]));
                    if(k==0){
                        P1=(PN[k]-PN[k+1])/2;
                        P2=(PN[k+1]-PN[k+2])/2;
                        P3=(PN[k+1]+PN[k+2])/2;                        
                            //**
                        V_TB=-M*(((Ut1[k][i+1][j]-Ut1[k][i][j])
                            +(Vt1[k][i][j+1]-Vt1[k][i][j]))*P2     
                            +((Ut1[k+1][i+1][j]-Ut1[k+1][i][j])
                            +(Vt1[k+1][i][j+1]-Vt1[k+1][i][j]))*P1)/(dex*(P1+P2));          
                        V_TB=roundValue(V_TB, 3);
                        double PA3= PN[k+1]-P3;
                        double PA2= PN[k]-PN[k+1];
                        double PA1=PN[k]-P3;                        
                        //**1
                        Wt[k][i][j]=(-V_TB*PA1+(Ws[i][j]*PA3/PA2-Wt[k+1][i][j]*PA2/PA3))/(PA3/PA2-PA2/PA3);
                        Wt[k][i][j]=roundValue(Wt[k][i][j], 3);
                        //**2
                        PHI[k][i][j]=(PHI[k+1][i][j]+R/Ro*TV*Math.pow((PN[k]/P0), RC)*(-PA2));
                        PHI[k][i][j]=roundValue(PHI[k][i][j], 3);
                    } else if (k==K-1){
                        double PA0= PN[k]-PN[k-1];
                        //***1
                        PHI[k][i][j]=(PHI[k-1][i][j]+R/Ro*TV*Math.pow(PN[k]/P0, RC)*PA0);
                        PHI[k][i][j]=roundValue(PHI[k][i][j], 3);
                        Wt[k][i][j]=0;
                    }else{
                        double PP1= PN[k-1]-PN[k];
                        double PP2= PN[k-1]-PN[k];
                        double PP3= PN[k]-PN[k+1];
                        PHI[k][i][j]=((R/Ro*TV*Math.pow(PN[k]/P0, RC)*PP1+PHI[k-1][i][j]*PP3/PP2
                            -PHI[k+1][i][j]*PP2/PP3))/(PP3/PP2-PP2/PP3);
                        PHI[k][i][j]=roundValue(PHI[k][i][j], 3);
                        V_TB=-M*((Ut1[k][i+1][j]-Ut1[k][i][j])
                            +(Vt1[k][i][j+1]-Vt1[k][i][j])
                            +(Ut1[k+1][i+1][j]-Ut1[k+1][i][j])
                            +(Vt1[k+1][i][j+1]-Vt1[k+1][i][j]))/(2*dex);
                        V_TB=roundValue(V_TB, 3);
                        Wt[k][i][j]=(-V_TB*PP1+Wt[k-1][i][j]*PP3/PP2-Wt[k+1][i][j]*PP2/PP1)
                            /(PP3/PP2- PP2/PP3);
                        Wt[k][i][j]=roundValue(Wt[k][i][j], 3);
                    }
                }
            }
        }
        P4=(PN[1]+PN[2])/2;
        for (int i = 1; i < X-1; i++) {
            for (int j = 1; j < Y-1; j++) {
                double VTS= M*((Ut1[1][i+1][j]-Ut1[1][i][j])+(Vt1[1][i][j+1]-Vt1[1][i][j]))/dex;
                Ws[i][j]=Ws[i][j]+VTS*(P4- PN[1]);
            }
        }
    }
    public static void printResult(){   
        System.out.println("U=");
        print3D(K, X, Y, Ut2);
        saveFile3D("Ut2.txt", Ut2);
        System.out.println("V=");
        saveFile3D("Vt2.txt", Vt2);
        print3D(K, X, Y, Vt2);
        System.out.println("TE=");
        saveFile3D("TE_t2.txt", TE_t2);
        print3D(K, X, Y, TE_t2);
        System.out.println("Q=");
        saveFile3D("Q_t2.txt", Q_t2);
        print3D(K, X, Y, Q_t2);
        System.out.println("Wt=");
        saveFile3D("Wt.txt", Wt);
        print3D(K, X, Y, Wt);
        System.out.println("PHI=");
        saveFile3D("PHI.txt", PHI);
        print3D(K, X, Y, PHI);
        
        System.out.println("PS=");
        print2D(X, Y, PS_t2);
        saveFile2D("PS_t2.txt", PS_t2);
        System.out.println("Ws=");
        saveFile2D("Ws.txt", Ws);
        print2D(X, Y, Ws);        
    }
    public static void saveFile3D(String fileName,double[][][] array ){
        try {            
            BufferedWriter bw= null;
            File file= new File(dataPath+fileName);
            if(!file.exists()){
                file.createNewFile();
            }
            bw= new BufferedWriter(new FileWriter(file));
            for (int i = 0; i < array.length; i++) {
                for (int j = 0; j < array[1].length; j++) {
                    for (int k = 0; k < array[0][1].length; k++) {
                        bw.write((i+1)+"\t"+(j+1)+"\t"+(k+1)+"\t"+array[i][j][k]);
                        bw.newLine();
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
        }
    }
     public static void saveFile2D(String fileName,double[][] array ){
        try {            
            BufferedWriter bw= null;
            File file= new File(dataPath+fileName);
            if(!file.exists()){
                file.createNewFile();
            }
            bw= new BufferedWriter(new FileWriter(file));
            for (int i = 0; i < array.length; i++) {
                for (int j = 0; j < array[1].length; j++) {                    
                        bw.write((i+1)+"\t"+(j+1)+"\t"+array[i][j]);
                        bw.newLine();                    
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
        }
    }
    public static void print3D(int k, int i, int j, double[][][] array){
        for (int l = 0; l < k; l++) {
            for (int m = 0; m < i; m++) {
                for (int n = 0; n < j; n++) {
                    System.out.println((l+1)+"\t"+(m+1)+"\t"+(n+1)+"\t"+array[l][m][n]);
                }
            }
        }
    }
    public static void print2D(int i, int j, double[][] array){
        for (int k = 0; k < i; k++) {
            for (int l = 0; l < j; l++) {
                System.out.println((k+1)+"\t"+(l+1)+"\t"+array[k][l]);
            }
        }
    }
    public static void saveResult(){
        
    }
    public static double roundValue(double value, int n){
        double decimal_place= Math.pow(10, n);
        return Math.round(value*decimal_place)/decimal_place;
    }
    public static void main(String[] args) {
        // TODO code application logic here
       
        System.out.println("CHUONG TRINH DU BAO VA CHUAN DOAN THOI TIET");
        System.out.println("-------------------------------------------");
        init();
        for(int n=0; n<N; n++){
            det1=2*det;
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
        printResult();
        
    }
}
