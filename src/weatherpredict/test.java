/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package weatherpredict;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 *
 * @author ngo
 */
public class test {
    public static void main(String[] args) throws FileNotFoundException {
        System.out.println(Math.pow(1000/25, 0.286));
        System.out.println(WeatherPredict.roundValue(5*Math.cos(Math.PI/12), 3));
    }    
}
