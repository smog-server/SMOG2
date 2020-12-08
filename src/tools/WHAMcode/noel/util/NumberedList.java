package noel.util;

import java.util.*;


public class NumberedList extends ArrayList implements Comparable<NumberedList> {
	public int num;
	public NumberedList(int num) {
		this.num = num;
	}
	public int compareTo(NumberedList list) {
	    if(num > list.num) return 1;
	    else if(num < list.num) return -1;
	    else return 0; //equal
	}
	public void print() {
	    Object[] lists = this.toArray();
	    for(int i = 0; i < lists.length; i++) {
	        System.out.print(((NumberedList)lists[i]).num+" ");
	    }
	    System.out.println();
	}
}   	