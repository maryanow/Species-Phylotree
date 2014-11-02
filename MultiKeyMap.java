import java.util.Map;
import java.lang.Double;

public class MultiKeyMap<V> {
    private java.util.HashMap<String,Double> map;

    public MultiKeyMap() {
        this.map = new java.util.HashMap<String,Double>();
    }

    public void put(String k1, String k2, double value) {
        if( k1.contains("|") || k2.contains("|") ) {
            System.err.println("Error: Keys in a MultiKeyMap cannot contain the bar character (\"|\")");
            System.exit(10);
        }
        map.put(k1 + "|" + k2, value);
        if( !k1.equals(k2) ) {
        	map.remove(k2 + "|" + k1);
        }
        return;
    }

    public Double get(String k1, String k2) {
        if( k1.contains("|") || k2.contains("|") ) {
            System.err.println("Error: Keys in a MultiKeyMap cannot contain the bar character (\"|\")");
            System.exit(10);
        }

        if (map.containsKey(k1 + "|" + k2)) {
            return map.get(k1 + "|" + k2);
        } 
	else {
            return map.get(k2 + "|" + k1);
        }
    }

    public void remove(String k1, String k2) {
        if (k1.contains("|") || k2.contains("|")) {
            System.err.println("Error: Keys in a MultiKeyMap cannot contain the bar character (\"|\")");
            System.exit(10);
        }

        map.remove(k1 + "|" + k2);
        map.remove(k2 + "|" + k1);
    }
  
    public String minDist(){
    	String s = "";
    	Double d = Double.POSITIVE_INFINITY;
    	for (Map.Entry<String, Double> entry : map.entrySet()){
   	    if (Double.compare(d, entry.getValue()) > 0){       // d is greater than entry.getValue()
    		d = entry.getValue();
    		s = entry.getKey();
    	    }
    	}

    	return s + "|" + Double.toString(d);
    }
}
