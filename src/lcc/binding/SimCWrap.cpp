/**
 *
 * C-wrapper for SimModule
 *
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>

class FSimModule {
public:
    FSimModule() 
        : SimModuleObject(nullptr) 
    {
        Py_Initialize();
        std::vector<std::string> Paths{"."};
        AddPath(Paths);

        SimModuleObject = LoadModule("SimModule");
    }

    ~FSimModule() {
        Py_XDECREF(SimModuleObject);
        Py_FinalizeEx();
    }

public:

    void Initiailze()
    {
        StateObject = CallFunction(SimModuleObject, "FState");
        DataSetObject = CallFunction(SimModuleObject, "FDataset");
        DataManagerObject = CallFunction(SimModuleObject, "FDataManager");

        {
            // Call SimModule.FSimulation(State, Data, DataManager)
            //
            PyObject *ArgsObject;
            ArgsObject = PyTuple_New(3);
            PyTuple_SetItem(ArgsObject, 0, StateObject);
            PyTuple_SetItem(ArgsObject, 1, DataSetObject);
            PyTuple_SetItem(ArgsObject, 2, DataManagerObject);

            SimObject = CallFunction(SimModuleObject, "FSimulation", ArgsObject);

            Py_DECREF(ArgsObject);
        }

        // Sim.Initialize()
        CallFunction(SimObject, "Initialize");
    }

    void RunOneStep()
    {
        assert(SimObject != NULL);

        CallFunction(SimObject, "SimLoop_WithoutSpatialSimulation");
        CallFunction(SimObject, "TriggerEventMoleculeCount");
        CallFunction(SimObject, "ExportData");
    }

    void GetLegend(std::vector<std::string>& OutLegend)
    {
        PyObject *LegendObject = PyObject_GetAttrString(DataSetObject, "Legend");

        if (LegendObject) {
            int LegendCount = PyList_Size(LegendObject);
            for(int i = 0; i < LegendCount; i++) {
                PyObject *p = PyList_GetItem(LegendObject, i);
                PyObject *s = PyUnicode_AsUTF8String(p);
                OutLegend.push_back(PyBytes_AsString(s));
                Py_XDECREF(s);
            }
        }

        Py_XDECREF(LegendObject);
    }

    void GetData(std::vector<float>& OutData)
    {
        PyObject *DataObject = CallFunction(DataSetObject, "ExportDataAsList");

        if (DataObject) {
            int DataCount = PyList_Size(DataObject);
            fprintf(stderr, "DataCount : %d\n", DataCount);
            for(int i = 0; i < DataCount; i++) {
                PyObject *p = PyList_GetItem(DataObject, i);
                OutData.push_back(PyFloat_AsDouble(p));
            }
        } else {

            std::cerr << "No Data object" << std::endl;

        }

        Py_XDECREF(DataObject);
    }

private:
    PyObject* SimModuleObject;

    PyObject* StateObject;
    PyObject* DataSetObject;
    PyObject* DataManagerObject;
    PyObject* SimObject;


private:
    PyObject* LoadModule(const std::string& InModuleName)
    {
        PyObject *pName = PyUnicode_DecodeFSDefault(InModuleName.c_str());
        if (pName == NULL) {
            std::cerr << "Can't load module: " << InModuleName << std::endl;
            return nullptr;
        }

        PyObject *pModule = PyImport_Import(pName);
        Py_DECREF(pName);

        return pModule;
    }


    void AddPath(const std::vector<std::string>& Paths)
    {
        PyObject *pSys = PyImport_ImportModule("sys");
        PyObject *pPath = PyObject_GetAttrString(pSys, "path");

        for(auto& path: Paths) {
            PyList_Append(pPath, PyUnicode_FromString(path.c_str()));
        }

        Py_DECREF(pPath);
        Py_DECREF(pSys);
    }

    PyObject* CallFunction(PyObject *InModuleObject, const char* InFuncName, PyObject *InArgsObject = NULL)
    {
        PyObject* FuncObject;
        PyObject* ReturnObject = NULL;

        assert(InModulePtr != NULL);

        FuncObject = PyObject_GetAttrString(InModuleObject, InFuncName);

        if (FuncObject && PyCallable_Check(FuncObject)) {

            PyObject* ArgsObject;

            if (InArgsObject == NULL) {
                ArgsObject = PyTuple_New(0);
            } else {
                ArgsObject = InArgsObject;
            }

            ReturnObject = PyObject_CallObject(FuncObject, ArgsObject);

            if (InArgsObject == NULL) {
                Py_DECREF(ArgsObject);
            }

            if (ReturnObject == NULL) {
                PyErr_Print();
                std::cerr << "Call failed \"" << InFuncName << "\"" << std::endl;
            }
        } else {
            if (PyErr_Occurred()) {
                PyErr_Print();
            }
            std::cerr << "Can't find function \"" << InFuncName << "\"" << std::endl;
        }

        Py_XDECREF(FuncObject);

        return ReturnObject;
    }
public:
};


int main(int argc, char *argv[])
{

    FSimModule Sim;

    Sim.Initiailze();
    Sim.RunOneStep();

    std::vector<std::string> Legend;
    Sim.GetLegend(Legend);
    for(auto& s: Legend) {
        fprintf(stderr, "%s\t", s.c_str());
    }
    fprintf(stderr, "\n");

    std::vector<float> Data;
    Sim.GetData(Data);
    for(auto& d: Data) {
        fprintf(stderr, "%f\t", d);
    }
    fprintf(stderr, "\n");

    return 0;
}
