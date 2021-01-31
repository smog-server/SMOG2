
package org.smogserver.io;

import java.io.*;
import java.util.*;
import org.smogserver.util.exception.*;
import java.lang.*;

/**
 * Has file reading/writing capabilities 
 * Created on August 7, 2003, 11:48 PM
 * @author  jknoel
 * @version 0.9
 */
public class FileIO extends java.io.File {

    private FileWriter writer;
    private FileReader reader;
    private BufferedReader bufReader;
    private StreamTokenizer tok;
    private int type;
	public boolean alreadyExisted = false;
    public final static int READING = 0;
    public final static int WRITING = 2;
    public final static int NO_OVER_WRITING = 4;
    public final static int BUFFERED_READING = 1;
    public final static int TOKEN_READING = 3;
    /**
     * Creates a READING FileIO object
     * @param filename file to read from
     */
    public FileIO(String filename) {
	    super(filename);
	    type = WRITING;
	    try {
		    writer = new FileWriter(this, false);
	    } catch (Exception e) { e.printStackTrace(); }
    }
    /**
     * Creates a FileIO object.  The types are as follows:
     * <li>FileIO.READING use of {@link #read(int numChar)}
     * <li>FileIO.BUFFERED_READING {@link #readLine()} and {@link #readAllLines(int skip)}
     * <li>FileIO.TOKEN_READING {@link #getTokens(Object[] tokens)}
     * <li>FileIO.WRITING {@link #write(String s)} 
	 * <li>FileIO.NO_OVER_WRITING {@link #write(String s)} (no FileWriter created if file already exists, sets alreadyExisted=true for error checking )
     */
    public FileIO(String filename, int type) {
	    super(filename);
	    if(!this.exists()&&type!=FileIO.WRITING&&type!=FileIO.NO_OVER_WRITING){
		    System.out.println("File: \""+filename+"\" does not exist!  Exiting.");
		    System.exit(1);
	    }    
	    this.type = type;
	    if (type == 0 || type == 1 || type == 3) {
		    try {
			    reader = new FileReader(this);
			    if (type == 1 || type == 3) {
				    bufReader = new BufferedReader(reader);
			    }
			    if (type == 3) {
				    tok = new StreamTokenizer(bufReader);
			    }
		    } catch (Exception e) { e.printStackTrace(); }
	    } else if (type == 2 || type == 4) { //makes a new file
			if(type == 4 && this.exists()) {
				this.alreadyExisted = true;
				writer = null;
			} else {
				//is there a parent directory?
				if(this.getParent() != null) {
					// if parent directory does not exist, make it
					if( ! (new File(this.getParent())).exists() ) {
						(new File(this.getParent())).mkdirs();
					}
				}
			    try {
				    writer = new FileWriter(this, false);
			    } catch (IOException e) { 
					System.out.println("Problem writing to file: "+filename);
					e.printStackTrace(); 
				}
			}
	    } else {
		    System.out.println("Bad argument sent to FileIO!");
	    }
    }
    /**
     * Changes a FileIO.WRITING to append.
     */
    public void appending() {
	    if (type == 2) {
		    try {
			    writer = new FileWriter(this, true);
		    } catch (Exception e) { e.printStackTrace(); }
	    }
    }
    /**
     * A FileIO.READING object reads a <code>numChar</code> number of characters.
     */
    public char[] read(int numChar) throws IOException, DoneException {
	    char[] chars = new char[numChar];
	    if (reader.read(chars) == -1) throw new DoneException();
	    return chars;
    }
    /**
     * Returns the next line in FileIO.BUFFERED_READING object.  
     * Returns <code>null</code> at end of file.
     * @return next line
     */
    public String readLine() {
	    String line = "Error";
	    try {
		    line = bufReader.readLine();
	    } catch (Exception e) { e.printStackTrace(); }
	    return line;
    }
    /**
     * Writes a string to a FileIO.WRITING object.
     * @param string to write
	 * @return true if successufl
     */  
    public boolean write(String s) {
		if(writer == null) {
			return false;
		}
	    try {
		    writer.write(s);
		    writer.flush();
	    } catch (Exception e) { 
			e.printStackTrace(); 
			return false;
		}
		return true;
    }
    /**
     * Closes any type of FileIO object.
     */
    public void close() {
	    if (type == 2) {
		    try {
			    writer.close();
		    } catch (Exception e) { e.printStackTrace(); }
	    } else if (type == 0) {
		    try {
			    reader.close();
		    } catch (Exception e) { e.printStackTrace(); }
	    } else if (type == 1) {
		    try {
			    bufReader.close();
		    } catch (Exception e) { e.printStackTrace(); }
	    }
    }
    /**
     * Reads all the lines in a FileIO.BUFFERED_READING object.  Skips designated number of lines.
     * @param skip number of lines to skip (at beginning)
     * @return lines in the file
     */    
    public String[] readAllLines(int skip) {
	    if (type == BUFFERED_READING) {
		    return readFromFile(this.getName(),skip);
	    } else {
		    return null;
	    }
    }  
    public static String[] readFromFile(String filename, int skip) {
	    String[] linesInFile = new String[10];
	    int totLines = 0, indexes = linesInFile.length;
	    File file = new File(filename);
	    String line = null;
	    if (file.exists())
	    { //make sure file exists
		    try
		    {
			    BufferedReader inFromFile = new BufferedReader(
					    new FileReader( file ) );

			    for (int i = -1; i < skip; i++) {
				    line = inFromFile.readLine();
			    }
			    int num = 0;
			    StringBuffer sB;
			    while (line != null)
			    {  //read in entire file
				    if (indexes - 1 == totLines)
				    {
					    linesInFile = enlarge(linesInFile);
					    indexes = linesInFile.length;
				    }
				    num = line.length();
				    sB = new StringBuffer(line);
				    linesInFile[totLines] = line;
				    totLines++;
				    line = inFromFile.readLine();
			    }
			    inFromFile.close(); //close the reader
		    }
		    catch (FileNotFoundException e)
		    {
			    System.out.println("File does not exist!");
		    }
		    catch (IOException e)
		    {
			    System.out.println("Problem with readLine method"+
					    "- last line read in-> "+line);
			    linesInFile = null;
			    return linesInFile;
		    }
		    catch (Exception e) { e.printStackTrace(); }
	    }
	    else {
		    System.out.println("File " + filename + " does not exist");
		    return linesInFile;
	    }
	    String[] temp = new String[totLines];
	    for (int j=0; j < totLines; j++) { temp[j] = linesInFile[j]; }
	    return temp; //return lines in the file
    }

    /**
     * Removes everything after a comment character in every line.  
     * Useful on the output of {@link #readAllLines(int skip)}.
     * @param lines of the file
     * @param commentChars the chars designating comments
     * @return edited lines
     */     
    public static String[] removeComments(String[] lines, char[] commentChars) {
	    int numLines = lines.length;
	    for (int i = 0; i < numLines; i++) {
		    String line = lines[i];
		    for (int j = 0; j < commentChars.length; j++) {
			    int index = line.indexOf(commentChars[j]);
			    if ( index > 0) line = line.substring(0,index-1);
			    if ( index == 0) line = "";
		    }
		    lines[i] = line;
	    }
	    return lines;
    }
    /**
     * Removes everything after a comment character in every line.  
     * Useful on the output of {@link #readAllLines(int skip)}.
     * @param lines of the file
     * @param commentChar the char designating comments
     * @return edited lines
     */     
    public static String[] removeComments(String[] lines, char commentChar) {
	    int numLines = lines.length;
	    for (int i = 0; i < numLines; i++) {
		    String line = lines[i];
		    int index = line.indexOf(commentChar);
		    if ( index > 0) line = line.substring(0,index-1);
		    if ( index == 0) line = "";
		    lines[i] = line;
	    }
	    return lines;
    }

    public static void main(String[] args) {
	    String[] lines = FileIO.readFromFile("C:\\Documents and Settings\\jknoel\\Desktop\\Results.txt",0);
	    char[] commentChars = {'!','$'};
	    lines = FileIO.removeComments(lines,commentChars);
	    lines = FileIO.getTokens(lines," ");
	    for (int i = 0; i < lines.length; i++) {
		    System.out.println(lines[i]);
	    }
    }
    //doubles array
    private static String[] enlarge(String[] f)
    {
	    int numLinesInArray = f.length;
	    String[] tempLines = new String[numLinesInArray * 2];
	    for (int i=0; i < numLinesInArray; i++)
	    {
		    tempLines[i] = f[i];
	    }
	    f = tempLines;
	    return f;
    }

    /**
     * Fills given array with next bunch of tokens and returns status.
     * status > 0 means that tokens[i > status] are bogus and EOF reached.
     * Only works for FileIO.TOKEN_READING.  If token is a number, parses 
     * as a double.
     * @param tokens array to fill with next tokens
     * @return int status
     */
    public int getTokens(Object[] tokens) throws IOException {
	    int count = 0, value = StreamTokenizer.TT_WORD, status = -1;
	    if (type == TOKEN_READING) {
		    while (tok.ttype != StreamTokenizer.TT_EOF && count < tokens.length) {
			    switch (tok.nextToken()) {
				    case StreamTokenizer.TT_NUMBER: tokens[count] = new Double(tok.nval);
								    break;
				    case StreamTokenizer.TT_WORD: tokens[count] = tok.sval;
								  break;
				    case StreamTokenizer.TT_EOF: status = count;
								 break;
			    }
			    count++;
		    }
	    }
	    return status;
    }
    /**
     * Breaks an array of strings into an array of tokens.
     * Maintains order.  Delimiter string " :" breaks tokens
     * on " " and ":".
     * @param file
     * @param delim delimiter string
     * @return array of tokens
     */ 
    public static String[] getTokens(String[] file, String delim) {
	    String[] allTokens = new String[5*file.length];
	    int lines = file.length;
	    int numTokens = 0;
	    for (int i=0; i < lines; i++) {
		    StringTokenizer tokens = new StringTokenizer(file[i], delim, false);
		    while (tokens.hasMoreTokens()) {
			    if (allTokens.length <= numTokens) { allTokens = enlarge(allTokens); }
			    allTokens[numTokens] = tokens.nextToken();
			    numTokens++;
		    }
	    }
	    String[] temp = new String[numTokens];
	    for (int j=0; j < numTokens; j++) { temp[j] = allTokens[j]; }
	    return temp;
    }
    /**
     * Breaks a string into an array of tokens.
     * Maintains order.  Delimiter string " :" breaks tokens
     * on " " and ":".
     * @param file
     * @param delim delimiter string
     * @return array of tokens
     */ 
    public static String[] getTokens(String file, String delim) {
	    String[] allTokens = new String[10];
	    int numTokens = 0;
	    StringTokenizer tokens = new StringTokenizer(file, delim, false);
	    while (tokens.hasMoreTokens()) {
		    if (allTokens.length <= numTokens) { allTokens = enlarge(allTokens); }
		    allTokens[numTokens] = tokens.nextToken();
		    numTokens++;
	    }
	    String[] temp = new String[numTokens];
	    for (int j=0; j < numTokens; j++) { temp[j] = allTokens[j]; }
	    return temp;
    }

    /**
     * Takes the output of getTokens and parses into Integer, if that fails
     * tries Double, if that fails keeps it a string.  Puts them in same
     * order into a vector.  Need to grab the values from the objects in your
     * code
     *@param s
     *@return tokens
     */
    public static Vector parseTokens(String[] s) {
	    Vector tokens = new Vector();
	    boolean notInt = false, notDouble = false;
	    for (int i = 0; i < s.length; i++) {
		    try {
			    tokens.add(Integer.valueOf(s[i]));
		    } catch (Exception e) { notInt = true; }
		    try {
			    if (notInt) tokens.add(Double.valueOf(s[i]));
		    } catch (Exception e) {notDouble = true; }
		    if (notDouble) tokens.add(s[i]);
		    notInt = notDouble = false;
	    }
	    return tokens;
    }
    /**
     * Prints the contents of a String[] to a file by specifing an
     * already created FileIO object
     *@param output
     *@param file
     */
    public static void printToFile(String[] output, FileIO file) {
	    for (int i=0; i < output.length;i++) {
		    if (output[i] == null) break;
		    file.write(output[i]);
	    }
    }
    /**
     * Prints the contents of a String[] to a file by specifing an
     * already created FileIO object
     *@param output
     *@param file
     */
    public static void printToFile(String output, FileIO file) {
	    if (file.type == READING) file.write(output);
	    else System.out.println("Not a READING FileIO object");
    }    
    public String toString() {
	    return super.toString();
    }
	
	/**
	* Wrapper for my use so I don't have to remember the right java.io commands
	*/
	public boolean isWriteable() {
		return true;
	}
		
}
