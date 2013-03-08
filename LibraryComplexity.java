/**
 * Taken from PicardTools 
 */
import java.text.NumberFormat;
class LibraryComplexity {

    public static void main(String[] args) {
	if (args.length != 2) {
	    System.out.println("Usage: java LibraryComplexity <reads> <unique reads>");
	    System.exit(0);
	}
	long readPairs = 0;
	long uniqueReadPairs = 0;
	try {
	    readPairs = Long.parseLong(args[0]);
	    uniqueReadPairs = Long.parseLong(args[1]);
	}
	catch (NumberFormatException e) {
	    System.err.println("Arguments must be integers");
	    System.exit(1);
	}
	long result;
	try {
			result = estimateLibrarySize(readPairs, uniqueReadPairs);
	}
	catch (NullPointerException e) {
			System.err.println("Library complexity undefined when total = " + readPairs + " and unique = " + uniqueReadPairs);
			return;
	}
	long result2 = readPairs*readPairs/(2*(readPairs-uniqueReadPairs));
	System.out.println("Library complexity (new): " + NumberFormat.getInstance().format(result));
	System.out.println("Library complexity (old): " + NumberFormat.getInstance().format(result2));
	
    }

/**
 * Estimates the size of a library based on the number of paired end molecules 
 * observed and the number of unique pairs observed.
 * <br>
 * Based on the Lander-Waterman equation that states:<br>
 *     C/X = 1 - exp( -N/X )<br>
 * where<br>
 *     X = number of distinct molecules in library<br>
 *     N = number of read pairs<br>
 *     C = number of distinct fragments observed in read pairs<br>
 */
    public static Long estimateLibrarySize(final long readPairs, 
					   final long uniqueReadPairs) {
	final long readPairDuplicates = readPairs - uniqueReadPairs;
	
	if (readPairs > 0 && readPairDuplicates > 0) {
	    long n = readPairs;
	    long c = uniqueReadPairs;
	    
	    double m = 1.0, M = 100.0;
	    
	    if (c >= n || f(m*c, c, n) < 0) {
		throw new IllegalStateException("Invalid values for pairs and unique pairs: " + n + ", " + c);
	    }

	    while( f(M*c, c, n) >= 0 ) {
		m = M;
		M *= 10.0;
	    }

	    double r = (m+M)/2.0;
	    double u = f( r * c, c, n );
	    int i = 0;

	    while (u != 0 && i < 1000) {
					if (u > 0) m = r;
					else M = r;
					r = (m+M)/2.0;
					u = f( r * c, c, n );
					i++;
	    }
	    if (i == 1000) {
					System.err.println("Iterated 1000 times, returning estimate");
	    }
	    
	    return (long) (c * (m+M)/2.0);
	}
	else {
	    return null;
	}
    }
    
    /** Method that is used in the computation of estimated library size. */
    private static double f(double x, double c, double n) {
	return c/x - 1 + Math.exp(-n/x);
    }
    
}    
