package main;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.UUID;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.Pair;

public class Molecule {
    String name;
    double x;
    double y;
    double z;
    //Temp variables used for checking if moves increase or decrease energy
    double tempx;
    double tempy;
    double tempz;
    double radius;
    private static final MersenneTwister r = new MersenneTwister();
    ArrayList<Atom> atoms;
    public Molecule(String n, double rad, ArrayList<Atom> a){
        name = n;
        radius = rad;
        atoms = a;
        resetTemps();
        setCenterOfMass(); //Also sets starting position, just in case of error
    }
    //For copying molecules from database List object
    public Molecule(Molecule m){
        name = m.name;
        radius = m.radius;
        x = m.x;
        y = m.y;
        z = m.z;
        atoms = new ArrayList<>();
        for (Atom a : m.atoms){
            atoms.add(new Atom(a));
        }
        resetTemps();
    }
    //Moves molecule to x,y,z position
    public void put(double x, double y, double z){
        this.x = x;
        this.y = y;
        this.z = z;
        for (Atom a : atoms){
            a.x = a.defX + x;
            a.y = a.defY + y;
            a.z = a.defZ + z;
        }
        resetTemps();
        setCenterOfMass();
    }
    //Rotates molecule randomly based on maxRot and saves values to temp vars to be used in energy calculations
    public void rotateTemp(double maxRot){
        if (atoms.size() <= 1) return;
        //Rotation amount as random double from -maxRot to maxRot
        double rotAmount = (r.nextDouble() - 0.5) * 2 * maxRot;
        int c = r.nextInt(3); //Choose rotation direction at random
        for (Atom a : atoms){ //Move atoms around center of molecule based on x, y, or z rotations
            a.x -= x;
            a.y -= y;
            a.z -= z;
            switch (c){
                case 0:
                    a.tempx = a.x;
                    a.tempy = Math.cos(rotAmount) * a.y + Math.sin(rotAmount) * a.z;
                    a.tempz = - 1* Math.sin(rotAmount) * a.y + Math.cos(rotAmount) * a.z;
                    break;
                case 1:
                    a.tempx = Math.cos(rotAmount) * a.x - Math.sin(rotAmount) * a.z;
                    a.tempy = a.y;
                    a.tempz = Math.sin(rotAmount) * a.x + Math.cos(rotAmount) * a.z;
                    break;
                case 2:
                    a.tempx = Math.cos(rotAmount) * a.x + Math.sin(rotAmount) * a.y;
                    a.tempy = -1 * Math.sin(rotAmount) * a.x + Math.cos(rotAmount) * a.y;
                    a.tempz = a.z;
                    break;
            }
            a.x += x;
            a.y += y;
            a.z += z;
            a.tempx += x;
            a.tempy += y;
            a.tempz += z;
        }
    }
    //Calculates energy for each atom in molecule based on full positions
    public double calcEnergy (Molecule m, HashMap<Pair<UUID, UUID>, double[]> pairwiseDbase) {
    	double energy = 0;
    	for (Atom a : atoms) {
    		for (Atom b : m.atoms) {
    			energy += a.calcEnergy(b, pairwiseDbase);
    		}
    	}
    	return energy;
    }
    //Calculates energy for each atom in molecule based on temp positions
    public double calcTempEnergy (Molecule m, HashMap<Pair<UUID, UUID>, double[]> pairwiseDbase) {
    	double energy = 0;
    	for (Atom a : atoms) {
    		for (Atom b : m.atoms) {
    			energy += a.calcTempEnergy(b, pairwiseDbase);
    		}
    	}
    	return energy;
    }
    //Called at end of any function where molecule is decided to not move or is moved with put(), sets temps to molecule's current position
    public void resetTemps() {
    	tempx = x;
    	tempy = y;
    	tempz = z;
    	for (Atom a : atoms) {
    		a.resetTemps();
    	}
    }
    //Moves molecule based on position of temps
    public void setTemps() {
    	x = tempx;
    	y = tempy;
    	z = tempz;
    	for (Atom a : atoms) {
    		a.setTemps();
    	}
    }
    //Calculates center of mass and centers molecule on it.
    public void setCenterOfMass(){
        double totMass = 0;
        double massPosX = 0;
        double massPosY = 0;
        double massPosZ = 0;
        for (Atom a : atoms){
            if (a.symbol.contains("*")){
                continue;
            }
            totMass += a.mass;
            massPosX += a.x * a.mass;
            massPosY += a.y * a.mass;
            massPosZ += a.z * a.mass;
        }
        this.x = massPosX / totMass;
        this.y = massPosY / totMass;
        this.z = massPosZ / totMass;

        resetTemps();
    }
}
