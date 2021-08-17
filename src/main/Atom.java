package main;
public class Atom {
	private static final double K_VALUE = 332.0587;
    String symbol;
    double x;
    double y;
    double z;
    //Def (default) variables used to save relative position to center of molecule; ensures that atoms are placed properly compared to center of molecule when propogating space
    double defX;
    double defY;
    double defZ;
    //Temp variables used for checking if moves increase or decrease energy
    double tempx;
    double tempy;
    double tempz;
    //a,b,c,d,q constants used for energy calculations
    double a;
    double b;
    double c;
    double d;
    double q;

    double mass;
    public Atom(String s, double x, double y, double z, double A, double B, double C, double D, double Q, double mass){
        symbol = s;
        this.x = x;
        this.y = y;
        this.z = z;
        defX = x;
        defY = y;
        defZ = z;
        resetTemps();
        a = A;
        b = B;
        c = C;
        d = D;
        q = Q;
        this.mass = mass;
    }
    //For copying molecules from database List object, called by Molecule.Molecule(Molecule m)
    public Atom (Atom atom){
        symbol = atom.symbol;
        x = atom.x;
        y = atom.y;
        z = atom.z;
        defX = atom.defX;
        defY = atom.defY;
        defZ = atom.defZ;
        resetTemps();
        a = atom.a;
        b = atom.b;
        c = atom.c;
        d = atom.d;
        q = atom.q;
        mass = atom.mass;
    }
    //Calculates energy based on atom full position
    public double calcEnergy(Atom atom) { //Calculates energy compared to rest of system, instead of total energy of system
    	double A = Math.sqrt(this.a * atom.a);
    	double B = (this.b + atom.b) / 2;
    	double C = Math.sqrt(this.c * atom.c);
    	double D = Math.sqrt(this.d * atom.d);
    	double r = Math.sqrt(Math.pow(atom.x - x, 2) + Math.pow(atom.y - y, 2) + Math.pow(atom.z - z, 2)); //Distance between atoms
    	return (K_VALUE * this.q * atom.q / r) + A * Math.exp(-1 * B * r) - (C / Math.pow(r,  6)) + (D / Math.pow(r,  12));
    }
   // Calculates energy based on atom temp position
    public double calcTempEnergy(Atom atom) { //Calculates energy compared to rest of system, instead of total energy of system
    	double A = Math.sqrt(this.a * atom.a);
    	double B = (this.b + atom.b) / 2;
    	double C = Math.sqrt(this.c * atom.c);
    	double D = Math.sqrt(this.d * atom.d);
    	double r = Math.sqrt(Math.pow(atom.x - tempx, 2) + Math.pow(atom.y - tempy, 2) + Math.pow(atom.z - tempz, 2)); //Distance between atoms
    	return (K_VALUE * this.q * atom.q / r) + A * Math.exp(-1 * B * r) - (C / Math.pow(r,  6)) + (D / Math.pow(r,  12));
    }
    //Sets temps to atom's current position, called by Molecule.resetTemps()
    public void resetTemps() {
    	tempx = x;
    	tempy = y;
    	tempz = z;
    }
    //Moves molecule based on position of temps, called by Molecule.setTemps()
    public void setTemps() {
    	x = tempx;
    	y = tempy;
    	z = tempz;
    }
}
