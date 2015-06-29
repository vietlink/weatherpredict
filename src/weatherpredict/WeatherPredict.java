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
    int[] PN= {1000, 900, 650, 400, 150};
    private static final int P0=1000;
    
    float[][][] U_t= new float[X][Y][K];
    float[][][] U_t1= new float[X][Y][K];
    float[][][] U_t2= new float[X][Y][K];
    
    float[][][] V_t= new float[X][Y][K];
    float[][][] V_t1= new float[X][Y][K];
    float[][][] V_t2= new float[X][Y][K];
    
    float[][][] TE= new float[X][Y][K];
    float[][][] TE1= new float[X][Y][K];
    float[][][] TE2= new float[X][Y][K];
    
    float[][][] Q_t= new float[X][Y][K];
    float[][][] Q_t1= new float[X][Y][K];
    float[][][] Q_t2= new float[X][Y][K];
    
    float[][][] PS= new float[X][Y][K];
    float[][][] PS1= new float[X][Y][K];
    float[][][] PS2= new float[X][Y][K];
    
    float[][][] Wt= new float[X][Y][K];
    float[][][] Ws= new float[X][Y][1];
    float[][][] PHI= new float[X][Y][K];
    
    public static void main(String[] args) {
        // TODO code application logic here
        
        for(int n=0; n<N; n++){
            for(int k=0; k<K; k++){
                
            }
        }
    }
    
    
}
