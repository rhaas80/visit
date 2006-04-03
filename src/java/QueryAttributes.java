package llnl.visit;

import java.util.Vector;
import java.lang.Double;
import java.lang.Integer;

// ****************************************************************************
// Class: QueryAttributes
//
// Purpose:
//    This class contains attributes used for query.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Fri Mar 17 15:08:37 PST 2006
//
// Modifications:
//   
// ****************************************************************************

public class QueryAttributes extends AttributeSubject
{
    // Enum values
    public final static int ELEMENTTYPE_ZONE = 0;
    public final static int ELEMENTTYPE_NODE = 1;

    public final static int VARTYPE_MESH = 0;
    public final static int VARTYPE_SCALAR = 1;
    public final static int VARTYPE_VECTOR = 2;
    public final static int VARTYPE_TENSOR = 3;
    public final static int VARTYPE_SYMMETRIC_TENSOR = 4;
    public final static int VARTYPE_ARRAY = 5;
    public final static int VARTYPE_LABEL = 6;
    public final static int VARTYPE_MATERIAL = 7;
    public final static int VARTYPE_SPECIES = 8;
    public final static int VARTYPE_CURVE = 9;
    public final static int VARTYPE_UNKNOWN = 10;

    public final static int DATATYPE_ACTUALDATA = 0;
    public final static int DATATYPE_ORIGINALDATA = 1;


    public QueryAttributes()
    {
        super(15);

        name = new String("");
        variables = new Vector();
        variables.addElement(new String("default"));
        resultsMessage = new String("");
        worldPoint = new double[3];
        worldPoint[0] = 0;
        worldPoint[1] = 0;
        worldPoint[2] = 0;
        domain = -1;
        element = -1;
        resultsValue = new Vector();
        resultsValue.addElement(new Double(0));
        elementType = ELEMENTTYPE_ZONE;
        timeStep = 0;
        varTypes = new Vector();
        dataType = DATATYPE_ACTUALDATA;
        pipeIndex = -1;
        useGlobalId = false;
        xUnits = new String("");
        yUnits = new String("");
    }

    public QueryAttributes(QueryAttributes obj)
    {
        super(15);

        int i;

        name = new String(obj.name);
        variables = new Vector(obj.variables.size());
        for(i = 0; i < obj.variables.size(); ++i)
            variables.addElement(new String((String)obj.variables.elementAt(i)));

        resultsMessage = new String(obj.resultsMessage);
        worldPoint = new double[3];
        worldPoint[0] = obj.worldPoint[0];
        worldPoint[1] = obj.worldPoint[1];
        worldPoint[2] = obj.worldPoint[2];

        domain = obj.domain;
        element = obj.element;
        resultsValue = new Vector(obj.resultsValue.size());
        for(i = 0; i < obj.resultsValue.size(); ++i)
        {
            Double dv = (Double)obj.resultsValue.elementAt(i);
            resultsValue.addElement(new Double(dv.doubleValue()));
        }

        elementType = obj.elementType;
        timeStep = obj.timeStep;
        varTypes = new Vector();
        for(i = 0; i < obj.varTypes.size(); ++i)
        {
            Integer iv = (Integer)obj.varTypes.elementAt(i);
            varTypes.addElement(new Integer(iv.intValue()));
        }
        dataType = obj.dataType;
        pipeIndex = obj.pipeIndex;
        useGlobalId = obj.useGlobalId;
        xUnits = new String(obj.xUnits);
        yUnits = new String(obj.yUnits);

        SelectAll();
    }

    public boolean equals(QueryAttributes obj)
    {
        int i;

        // Compare the worldPoint arrays.
        boolean worldPoint_equal = true;
        for(i = 0; i < 3 && worldPoint_equal; ++i)
            worldPoint_equal = (worldPoint[i] == obj.worldPoint[i]);

        // Create the return value
        return ((name == obj.name) &&
                (variables == obj.variables) &&
                (resultsMessage == obj.resultsMessage) &&
                worldPoint_equal &&
                (domain == obj.domain) &&
                (element == obj.element) &&
                (resultsValue == obj.resultsValue) &&
                (elementType == obj.elementType) &&
                (timeStep == obj.timeStep) &&
                (varTypes == obj.varTypes) &&
                (dataType == obj.dataType) &&
                (pipeIndex == obj.pipeIndex) &&
                (useGlobalId == obj.useGlobalId) &&
                (xUnits == obj.xUnits) &&
                (yUnits == obj.yUnits));
    }

    // Property setting methods
    public void SetName(String name_)
    {
        name = name_;
        Select(0);
    }

    public void SetVariables(Vector variables_)
    {
        variables = variables_;
        Select(1);
    }

    public void SetResultsMessage(String resultsMessage_)
    {
        resultsMessage = resultsMessage_;
        Select(2);
    }

    public void SetWorldPoint(double[] worldPoint_)
    {
        worldPoint[0] = worldPoint_[0];
        worldPoint[1] = worldPoint_[1];
        worldPoint[2] = worldPoint_[2];
        Select(3);
    }

    public void SetWorldPoint(double e0, double e1, double e2)
    {
        worldPoint[0] = e0;
        worldPoint[1] = e1;
        worldPoint[2] = e2;
        Select(3);
    }

    public void SetDomain(int domain_)
    {
        domain = domain_;
        Select(4);
    }

    public void SetElement(int element_)
    {
        element = element_;
        Select(5);
    }

    public void SetResultsValue(Vector resultsValue_)
    {
        resultsValue = resultsValue_;
        Select(6);
    }

    public void SetElementType(int elementType_)
    {
        elementType = elementType_;
        Select(7);
    }

    public void SetTimeStep(int timeStep_)
    {
        timeStep = timeStep_;
        Select(8);
    }

    public void SetVarTypes(Vector varTypes_)
    {
        varTypes = varTypes_;
        Select(9);
    }

    public void SetDataType(int dataType_)
    {
        dataType = dataType_;
        Select(10);
    }

    public void SetPipeIndex(int pipeIndex_)
    {
        pipeIndex = pipeIndex_;
        Select(11);
    }

    public void SetUseGlobalId(boolean useGlobalId_)
    {
        useGlobalId = useGlobalId_;
        Select(12);
    }

    public void SetXUnits(String xUnits_)
    {
        xUnits = xUnits_;
        Select(13);
    }

    public void SetYUnits(String yUnits_)
    {
        yUnits = yUnits_;
        Select(14);
    }

    // Property getting methods
    public String   GetName() { return name; }
    public Vector   GetVariables() { return variables; }
    public String   GetResultsMessage() { return resultsMessage; }
    public double[] GetWorldPoint() { return worldPoint; }
    public int      GetDomain() { return domain; }
    public int      GetElement() { return element; }
    public Vector   GetResultsValue() { return resultsValue; }
    public int      GetElementType() { return elementType; }
    public int      GetTimeStep() { return timeStep; }
    public Vector   GetVarTypes() { return varTypes; }
    public int      GetDataType() { return dataType; }
    public int      GetPipeIndex() { return pipeIndex; }
    public boolean  GetUseGlobalId() { return useGlobalId; }
    public String   GetXUnits() { return xUnits; }
    public String   GetYUnits() { return yUnits; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(name);
        if(WriteSelect(1, buf))
            buf.WriteStringVector(variables);
        if(WriteSelect(2, buf))
            buf.WriteString(resultsMessage);
        if(WriteSelect(3, buf))
            buf.WriteDoubleArray(worldPoint);
        if(WriteSelect(4, buf))
            buf.WriteInt(domain);
        if(WriteSelect(5, buf))
            buf.WriteInt(element);
        if(WriteSelect(6, buf))
            buf.WriteDoubleVector(resultsValue);
        if(WriteSelect(7, buf))
            buf.WriteInt(elementType);
        if(WriteSelect(8, buf))
            buf.WriteInt(timeStep);
        if(WriteSelect(9, buf))
            buf.WriteIntVector(varTypes);
        if(WriteSelect(10, buf))
            buf.WriteInt(dataType);
        if(WriteSelect(11, buf))
            buf.WriteInt(pipeIndex);
        if(WriteSelect(12, buf))
            buf.WriteBool(useGlobalId);
        if(WriteSelect(13, buf))
            buf.WriteString(xUnits);
        if(WriteSelect(14, buf))
            buf.WriteString(yUnits);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetName(buf.ReadString());
                break;
            case 1:
                SetVariables(buf.ReadStringVector());
                break;
            case 2:
                SetResultsMessage(buf.ReadString());
                break;
            case 3:
                SetWorldPoint(buf.ReadDoubleArray());
                break;
            case 4:
                SetDomain(buf.ReadInt());
                break;
            case 5:
                SetElement(buf.ReadInt());
                break;
            case 6:
                SetResultsValue(buf.ReadDoubleVector());
                break;
            case 7:
                SetElementType(buf.ReadInt());
                break;
            case 8:
                SetTimeStep(buf.ReadInt());
                break;
            case 9:
                SetVarTypes(buf.ReadIntVector());
                break;
            case 10:
                SetDataType(buf.ReadInt());
                break;
            case 11:
                SetPipeIndex(buf.ReadInt());
                break;
            case 12:
                SetUseGlobalId(buf.ReadBool());
                break;
            case 13:
                SetXUnits(buf.ReadString());
                break;
            case 14:
                SetYUnits(buf.ReadString());
                break;
            }
        }
    }


    // Attributes
    private String   name;
    private Vector   variables; // vector of String objects
    private String   resultsMessage;
    private double[] worldPoint;
    private int      domain;
    private int      element;
    private Vector   resultsValue; // vector of Double objects
    private int      elementType;
    private int      timeStep;
    private Vector   varTypes; // vector of Integer objects
    private int      dataType;
    private int      pipeIndex;
    private boolean  useGlobalId;
    private String   xUnits;
    private String   yUnits;
}

