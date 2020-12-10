package org.smogserver.scm;

import org.smogserver.io.*;

//arg[0] = gromacs bonds vs time distance.xvg file
//arg[1] = native distance.xvg file
//arg[2] = output cont
//arg[3] = output cont i
//arg[4] = the gamma!
public class Contacts {

	public static void main(String[] args){
		FileIO file1 = new FileIO(args[0],FileIO.BUFFERED_READING);
		FileIO file2 = new FileIO(args[1],FileIO.BUFFERED_READING);
		String line = file2.readLine();
		String[] dist = FileIO.getTokens(line," ");
		int numC = dist.length-1;
		double[] nat = new double[numC];
		double gamma = Double.parseDouble(args[4]);
		for(int i = 1; i < numC+1; i++){
			nat[i-1]=Double.parseDouble(dist[i])*gamma;
		}
		file2.close();

		FileIO out = new FileIO(args[2],FileIO.WRITING);
		FileIO outi = new FileIO(args[3],FileIO.WRITING);

		line=file1.readLine();
		while(line != null){
			dist = FileIO.getTokens(line," ");
			double d;
			int Q = 0;
			for(int i = 1; i < numC+1; i++){
				d=Double.parseDouble(dist[i]);
				if(d < nat[i-1]) {
					Q++;
					outi.write(i+"\n");
				}
			}
			out.write(Q+"\n");
			Q=0;
			line=file1.readLine();
		}
		file1.close();
		out.close();
		outi.close();
	}
}

