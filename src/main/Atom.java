package main;

import org.apache.commons.math3.util.Pair;
import java.util.HashMap;
import java.util.UUID;

public class Atom {
    public UUID uuid;

	public static final double K_VALUE = 332.0637133;
    String symbol;
    double x;
    double y;
    double z;
    //Def (default) variables used to save relative position to center of molecule; ensures that atoms are placed properly compared to center of molecule when propagating space
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
    //Composite of K_VALUE and q, to remove one floating point operation per energy calculation
    double kTimesQ;

    double mass;

    public Atom(String s, double x, double y, double z, double A, double B, double C, double D, double Q, double mass){
        this.uuid = UUID.randomUUID();
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
        this.kTimesQ = K_VALUE * q;
    }

    public Atom(String s, double... params){
        this.uuid = UUID.randomUUID();
        symbol = s;
        this.x = params[0];
        this.y = params[1];
        this.z = params[2];
        defX = x;
        defY = y;
        defZ = z;
        resetTemps();
        a = params[3];
        b = params[4];
        c = params[5];
        d = params[6];
        q = params[7];
        this.mass = params[8];
        this.kTimesQ = K_VALUE * q;
    }
    //For copying molecules from database List object, called by Molecule.Molecule(Molecule m)
    public Atom (Atom atom){
        this.uuid = atom.uuid;
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
        kTimesQ = K_VALUE * q;
    }
    //Calculates energy based on atom full position
    public double calcEnergy(Atom atom, HashMap<Pair<UUID, UUID>, double[]> pairwiseDbase) { //Calculates energy compared to rest of system, instead of total energy of system
        double A;
        double B;
        double C;
        double D;
        Pair<UUID, UUID> key = new Pair<>(this.uuid, atom.uuid);
        double[] value = pairwiseDbase.get(key);
        A = value[0];
        B = value[1];
        C = value[2];
        D = value[3];
    	double r = Math.sqrt(Math.pow(atom.x - x, 2) + Math.pow(atom.y - y, 2) + Math.pow(atom.z - z, 2)); //Distance between atoms
    	return (kTimesQ * atom.q / r) + A * Math.exp(-1 * B * r) - (C / Math.pow(r,  6)) + (D / Math.pow(r,  12));
    }
   // Calculates energy based on atom temp position
    public double calcTempEnergy(Atom atom, HashMap<Pair<UUID, UUID>, double[]> pairwiseDbase) { //Calculates energy compared to rest of system, instead of total energy of system
        double A;
        double B;
        double C;
        double D;
        Pair<UUID, UUID> key = new Pair<>(this.uuid, atom.uuid);
        double[] value = pairwiseDbase.get(key);
        A = value[0];
        B = value[1];
        C = value[2];
        D = value[3];
    	double r = Math.sqrt(Math.pow(atom.x - tempx, 2) + Math.pow(atom.y - tempy, 2) + Math.pow(atom.z - tempz, 2)); //Distance between atoms
    	return (kTimesQ * atom.q / r) + A * Math.exp(-1 * B * r) - (C / Math.pow(r,  6)) + (D / Math.pow(r,  12));
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
    @Override
    public final boolean equals(Object object){
        return object.getClass().equals(this.getClass()) && this.symbol.equals(((Atom) object).symbol);
    }
    @Override
    public final int hashCode(){
        int result = 17;
        result = result * 31 + symbol.hashCode();
        return result;
    }
}
