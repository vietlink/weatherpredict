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
        File f= new File("H:\\git\\PythonTut\\PythonTut\\src\\test\\test1.txt");
        Scanner s= new Scanner(f);
        int[] a= new int[5];
        for (int i = 0; i < a.length; i++) {
            a[i]=s.nextInt();
        }
        for (int b : a) {
            System.out.println(b);
        }
    }
}
