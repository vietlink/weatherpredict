/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package weatherpredict;

import java.io.File;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import java.util.Scanner;

/**
 *
 * @author ngo
 */
public class test {
    public static void main(String[] args) throws FileNotFoundException {
        DecimalFormat df= new DecimalFormat("#0.###");
        float  a= 1.000000443245f;
//        a= WeatherPredict.roundValue((float) a, 6);
        df.format(a);
        System.out.println(a);
    }    
}
