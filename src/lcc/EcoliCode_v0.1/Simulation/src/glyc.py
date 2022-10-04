import numpy as np
import warnings
from random import random, randrange, randint
from dataclasses import dataclass

from util import printTimestamp, printBlockMessage
from units import cnt2mol, mol2cnt
from modelEq import saturation, alloAct, alloInhib

EXCLUDE_FROM_PLOT = []#['HK', 'PEPWELD', 'PK']
ONE_COUNT_IN_MICROMOLAR = 1.66e-18

DICT_TIME = {
    'simDurationSeconds':5000, #172800, #48 hrs if 1s
    'simStepsPerSecond': 1,
    }
# Format: moleculename: (concentration uM, type) 
# Type can be 'met' for metabolite or 'enz' for enzyme
DICT_T0_CONC = {
            # Changeable met concs
            'Gluc' : (5000, 'met'),
            'F6P' :  (0, 'met'),          # Goal: 8.8mM - https://www.nature.com/articles/nchembio.186.pdf,
            'F16BP' :(0, 'met'),        # Goal: 15.0mM - https://www.nature.com/articles/nchembio.186.pdf
            'PEP' :  (0, 'met'),         # Goal: 0.18mM - https://www.nature.com/articles/nchembio.186.pdf
            'ADP' :  (560, 'met'),           # Goal: 0.56mM - https://www.nature.com/articles/nchembio.186.pdf
            'ATP' :  (9600, 'met'),
            'PYR' :  (0, 'met'),
            # Changeable enz concs
            'PFKa' : (1e-3, 'enz'),
            'PFKi' : (1e-3, 'enz'),
            'PFKb' : (5e-3, 'enz'),
            # Unchanging Enzyme Concs:
            'HK' :   (7e-3, 'enz'),
            'PEPWELD' : (4.5e-3, 'enz'),
            'PK' :   (7.0e-2, 'enz'),
        }
DICT_KCAT = {
    'PFKa':1e2,
    'PFKi':1e2,
    'HK' : 4e3,
    'PK' : 1.7e3,
    'PEPWELD' : 1e3
    }
DICT_K = {
    'PFKa-->PFKb':1e-3,
    'PFKi-->PFKb':1e-3,
    'PFKa-->PFKi':cnt2mol(5e11, 'micro'),
    'PFKi-->PFKa':cnt2mol(5e11, 'micro'),
    }
DICT_KM = {
    'PFKb_PFKb':5e-7,
    'PFKb_PFKa':1e-2,
    'PFKb_PFKi':1e-2,
    'PFKb_ADP': 550.0,
    'PFKb_PEP': 5e-2,
    'PFKa_F6P': 1.0,
    'PFKa_ATP': 60,
    'PFKi_F6P': 1e7,
    'PFKi_ATP': 60,
    'HK_Gluc' : 500.0,
    'HK_ATP' : 15000.0,
    'PEPWELD_F16BP' : 5000.0,
    'PK_PEP' : 35.0,
    'PK_ADP' : 50.0,
    }
DICT_K_REGULATORY = {
    # Format: Ka/Ki_Molecule_Target
    # Ki
    'Ki_PYR_HK' : 50.0,
    'Ki_PEP_PFK' : 100.0,
    'Ki_ATP_PK' : 12000.0,
    'Ki_PEP_PEPWELD' : 400.0,
    'Ki_F6P_HK' : 20.0,
    # Ka
    'Ka_F16BP_PK' : 100.0
}

DICT_ENZ_SUBSTRATE = {
    'HK' : {'subs':[('Gluc','sat'), ('ATP','sat')], 'prod':[('F6P','sat'), ('ADP','sat')], 'alloI':[('F6P','AlloI')],},
    'PFKa' :  {'subs':[('F6P','sat'), ('ATP','sat')], 'prod':[('F16BP','sat'), ('ADP','sat')]},
    'PFKi' :{'subs':[('F6P','sat'), ('ATP','sat')], 'prod':[('F16BP','sat'), ('ADP','sat')]},
    'PFKb' : {'subs':[('ADP','sat'), ('PEP','sat')], 'prod':[('PFKa','sat'), ('PFKi','sat')]},
    'PEPWELD' : {'subs':[('F16BP','sat')], 'prod':[('PEP','sat')]},
    'PK' : {'subs':[('PEP','sat'), ('ADP','sat')], 'prod':[('PYR','sat'), ('ATP','sat')], 'alloI':[('ATP', 'AlloI')], 'alloA':[('F16BP', 'AlloA')]},    
}

# Breaking this down into molecules, stochiometry, and reactions
Molecules = np.empty(shape = len(DICT_T0_CONC.keys()))




class Glycolysis:
    """
    TODO:
    -----
    Conserved quantities checks
    
    Notes:
    ------
     --> current task --> Current model does not have Glucose saturation in the PFK, just assumes it is always at saturation.
        
    """
    def __init__(self, dictTime, dictConc, dictKcat, dictK, dictKM, dictKReg, verbose:bool = True, debug:bool = True):
        # Time parameters
        self.time_simDurationSeconds = dictTime['simDurationSeconds']
        self.time_simStepsPerSecond = dictTime['simStepsPerSecond']
        self.time_totalSteps = self.time_simDurationSeconds * self.time_simStepsPerSecond
        self.time_stepResolution = 1 / self.time_simStepsPerSecond
        self.percComplete = 0

        # Adding in reaction params
        self.dictKcat = dictKcat
        self.dictK = dictK
        self.dictKM = dictKM
        self.dictKreg = dictKReg

        # Init Enzymes, metabolites, etc.
        self.concLegend = [key for key in dictConc.keys()]
        self.enzymeLegend = [key for key in dictConc.keys() if dictConc[key][1] == 'enz']
        # Init Containers for simulation values
        self.time = np.linspace(0, self.time_totalSteps, self.time_totalSteps)
        self.dictCountArrays= {i:np.ones([self.time_totalSteps]) for i in self.concLegend} 
        self.zerosCheck = {i:True for i in self.concLegend} 
        self.molecules = {}  

        self.enzymeMetaboliteLegend = [""]
        #self.enzymeMetaboliteLegend = [{'enz':key, 'subs':val} for key in DICT_ENZ_SUBSTRATE.keys() for val in DICT_ENZ_SUBSTRATE[key]]
        #print(self.enzymeMetaboliteLegend)

        # TODO: To be reimplemented
        self.dictEnzymePercentSubstrateSaturation = {i:np.zeros([self.time_totalSteps]) for i in self.enzymeMetaboliteLegend}
        self.dictEnzymePercentMaxActivty = {i:np.zeros([self.time_totalSteps]) for i in self.enzymeLegend} 
        self.dictRawStepTurnover = {i:np.zeros([self.time_totalSteps]) for i in self.concLegend}
        self.dictAdjustedStepTurnover = {i:np.zeros([self.time_totalSteps]) for i in self.concLegend}  
        self.dictRawEnzymeStepTurnover = {i:np.zeros([self.time_totalSteps]) for i in self.enzymeLegend}
        self.dictAdjustedEnzymeStepTurnover = {i:np.zeros([self.time_totalSteps]) for i in self.enzymeLegend}  
    
        # QSSA
        self.QSSAFlag = {f"{i}_{j[0]}":False for i in DICT_ENZ_SUBSTRATE.keys() for j in DICT_ENZ_SUBSTRATE[i]['subs']}
        self.QSSAThreshold = 0.01
        self.dictQSSACheck = {f"{i}_{j[0]}":np.zeros([self.time_totalSteps]) for i in DICT_ENZ_SUBSTRATE.keys() for j in DICT_ENZ_SUBSTRATE[i]['subs']}  

        # Init t0 Counts
        for key in self.concLegend:
            # The container for all molecule counts
            self.dictCountArrays[key][0] = dictConc[key][0]  
            # The container for current molecule counts
            self.molecules[key] = dictConc[key][0]

        self.debug = debug
        self.verbose = verbose

        self.passPYRConsumption = randint(1, 5) * 10 ** randrange(-3, -1, 1)
        self.passATPConsumption = randint(7, 8) * 10 ** randrange(-3, -1, 1)            
        self.passPEPConsumption = randint(1, 3) * 10 ** randrange(-3, -1, 1)  

    # Reaction Equations
    # Gluc --> F6P
    def hexokinase(self):
        """ Gluc + ATP --> F6P + ADP; Catalyzed by HK; Allosterically inhibited by Pyruvate and F6P"""
        vmax = self.dictKcat['HK'] * self.molecules['HK'] 
        sat = saturation(self.molecules['Gluc'], self.dictKM['HK_Gluc']) * saturation(self.molecules['ATP'], self.dictKM['HK_ATP'])
        allo =  alloInhib(self.molecules['F6P'], self.dictKreg['Ki_F6P_HK'])
        return vmax * sat * allo

    # PFKi <--> PFKb <--> PFKa
    def pfkActiveToBase(self):
        return self.dictK['PFKa-->PFKb'] * saturation(self.molecules['PFKa'], self.dictKM['PFKb_PFKa']) 
    def pfkInactiveToBase(self):
        return  self.dictK['PFKi-->PFKb'] * saturation(self.molecules['PFKi'], self.dictKM['PFKb_PFKi']) 
    def reversible_PFK(self):
        """ Reversible PFK (PFKa <--> PFKi) with Michaelis Menten Kinetics 1 Substrate 1 Product """
        num =  self.molecules['PFKb'] / self.dictKM['PFKb_PFKb'] * (
        self.dictK['PFKa-->PFKi'] * saturation(self.molecules['PFKa'], self.dictKM['PFKb_PFKa']) * saturation(self.molecules['ADP'],self.dictKM['PFKb_ADP']) -
        self.dictK['PFKi-->PFKa'] * saturation(self.molecules['PFKi'], self.dictKM['PFKb_PFKi']) * saturation(self.molecules['PEP'],self.dictKM['PFKb_PEP'])
        )
        denom = 1 + (self.molecules['PFKa'] / self.dictKM['PFKb_PFKa'] ) + (self.molecules['PFKi'] / self.dictKM['PFKb_PFKi'])
        return num / denom
    # F6P --> F16BP
    def pfk_active(self):
        """ PFKa Catalysis (F6P --> F16BP) Michaelis Menten 1 substrate irreversible"""
        vmax = self.dictKcat['PFKa'] * self.molecules['PFKa'] 
        sat = saturation(self.molecules['F6P'], self.dictKM['PFKa_F6P']) * saturation(self.molecules['ATP'], self.dictKM['PFKa_ATP'])
        allo = 1
        return vmax * sat * allo
    def pfk_inactive(self):
        """ PFKi Catalysis (F6P --> F16BP) Michaelis Menten 1 substrate irreversible"""
        vmax = self.dictKcat['PFKi'] * self.molecules['PFKi'] 
        sat = saturation(self.molecules['F6P'], self.dictKM['PFKi_F6P']) * saturation(self.molecules['ATP'], self.dictKM['PFKi_ATP'])
        allo = 1
        return vmax * sat * allo

    # F16BP --> PEP
    def pepweldase(self):
        vmax = self.dictKcat['PEPWELD'] * self.molecules['PEPWELD']
        sat = saturation(self.molecules['F16BP'], self.dictKM['PEPWELD_F16BP']) 
        allo = 1
        return vmax * sat * allo

    # PEP --> PYR
    def pyruvateKinase(self):
        vmax = self.dictKcat['PK'] * self.molecules['PK']
        sat = saturation(self.molecules['PEP'], self.dictKM['PK_PEP']) * saturation(self.molecules['ADP'], self.dictKM['PK_ADP'])
        allo = alloInhib(self.molecules['ATP'], self.dictKreg['Ki_ATP_PK']) * alloAct(self.molecules['F16BP'], self.dictKreg['Ka_F16BP_PK'])
        return vmax * sat * allo

    # Passive ATP Consumption:
    def passiveATPConsumption(self):
        return self.passATPConsumption * self.molecules['ATP']
    def passivePYRConsumption(self):
        return self.passPYRConsumption * self.molecules['PYR']
    def passivePEPConsumption(self):
        return self.passPEPConsumption * self.molecules['PEP']   

    # Concentration Change equations
    def dF6P_dt(self):
        influx = self.hexokinase()
        efflux = self.pfk_active() + self.pfk_inactive()
        return influx - efflux
    def dF16BP_dt(self):
        influx = self.pfk_active() + self.pfk_inactive()
        efflux = self.pepweldase()
        return influx - efflux
    def dPEP_dt(self):
        influx = self.pfkInactiveToBase() + self.pepweldase()
        efflux = abs(min(0, self.reversible_PFK())) + self.passivePEPConsumption()
        return influx - efflux
    
    def dATP_dt(self):
        influx = self.pyruvateKinase()
        efflux = self.hexokinase() + self.pfk_active() + self.pfk_inactive() + self.passiveATPConsumption()
        return influx - efflux 
    def dADP_dt(self):
        influx = self.pfkActiveToBase() + self.hexokinase() + self.pfk_active() + self.pfk_inactive() + self.passiveATPConsumption()
        efflux = max(0, self.reversible_PFK()) + self.pyruvateKinase()
        return influx - efflux
    
    def dPFKa_dt(self):
        influx = max(0, self.reversible_PFK())
        efflux = self.pfkActiveToBase()
        return influx - efflux
    def dPFKi_dt(self):
        influx = abs(min(0, self.reversible_PFK()))
        efflux = self.pfkInactiveToBase()
        return influx - efflux
    def dPFKb_dt(self):
        influx = self.pfkActiveToBase() + self.pfkInactiveToBase()
        efflux = abs(self.reversible_PFK()) # Basically always losing 
        return influx - efflux

    # 0 change values
    def dGluc_dt(self):
        influx = 0
        efflux = self.hexokinase() 
        return influx - efflux
    def dPYR_dt(self):
        influx = self.pyruvateKinase()
        efflux = self.passivePYRConsumption()
        return influx - efflux
    def dZeroChange_dt(self):
        return 0
    def rate(self):
        return {
            'Gluc' : self.dGluc_dt() * self.time_stepResolution,
            'F6P' : self.dF6P_dt() * self.time_stepResolution, 
            'F16BP' : self.dF16BP_dt() * self.time_stepResolution, 
            'PEP' : self.dPEP_dt() * self.time_stepResolution, 
            'ADP' : self.dADP_dt() * self.time_stepResolution,  
            'ATP' : self.dATP_dt() * self.time_stepResolution,  
            'PFKa' : self.dPFKa_dt() * self.time_stepResolution,
            'PFKi' : self.dPFKi_dt() * self.time_stepResolution,
            'PFKb' : self.dPFKb_dt() * self.time_stepResolution,
            'PYR' : self.dPYR_dt() * self.time_stepResolution,
            # Enzymes with fixed concentrations
            'HK' : self.dZeroChange_dt(),
            'PEPWELD' : self.dZeroChange_dt(),
            'PK' : self.dZeroChange_dt()
        }

    def debugSetRatesToZero(self, lstMolecules:list, dictRates:dict) -> dict:
        """Set specific molecule turnover rates to zero (for debugging)

        Parameters:
        -----------
            lstMoleclues : `list`
                List of molecule names which should have their rates set to 0. These values must match the keys of dictRates.
            dictRates : `dict` 
                Dictionary containing the step turn over (i.e. rate) and molecule name. format: {moleculeName : numberOfMoleculesToTurnOver} 
        Returns:
        --------
            zeroedRates : `dict`
                The input dictionary with zeroed rates matching the lstMolecules.
        """ 
        zeroedDictRates = {}
        for key in dictRates.keys():
            if key in lstMolecules:
                zeroedDictRates[key] = 0
            else:
                zeroedDictRates[key] = dictRates[key]
        return zeroedDictRates

    def updateMolecules(self, step):
        for key in self.molecules.keys():
            self.molecules[key] = self.dictCountArrays[key][step]

    def perturbation(self, step):
        # if step:
        #     newGluc = 5*random() * 10 ** randrange(-12, -9, 1)
        #     self.dictCountArrays['Gluc'][step] = newGluc
        #     self.molecules['Gluc'] = newGluc           
        
        #if step % 100 == 0:
            #self.passPYRConsumption = randint(1, 5) * 10 ** randrange(-2, -1, 1)
            #self.passATPConsumption = (randint(1,12) - 1/randint(3, 8)) * 10 ** randrange(-6, -5, 1)            
            #self.passPEPConsumption = randint(1, 3) * 10 ** randrange(-2, -1, 1)   

        if step % 500 == 0:
            newGluc = 5*random() * 10 ** randrange(0,6, 1)
            self.dictCountArrays['Gluc'][step] = newGluc
            self.molecules['Gluc'] = newGluc   

        pass

    def saveStepVmax(self, step):
        """ Look at each enzyme, get the current step vmax."""
        for key in self.dictEnzymePercentMaxActivty.keys():
            self.dictEnzymePercentMaxActivty[key][step] = self.dictKcat[key] * self.molecules[key]

    def checkAssertions(self,step):
        # Total mass conservation:
        totalPFK = self.dictCountArrays['PFKa'][0] + self.dictCountArrays['PFKi'][0] + self.dictCountArrays['PFKb'][0]
        if abs(self.molecules['PFKa'] + self.molecules['PFKi'] + self.molecules['PFKb'] - totalPFK) > 1e-12:
            if self.debug:
                warnings.warn(f"On Step: {step} PFK mass Not Conserved! \n Correct Total: {totalPFK}  Actual Total: {self.molecules['PFKa'] + self.molecules['PFKi'] + self.molecules['PFKb']}  Difference: {self.molecules['PFKa'] + self.molecules['PFKi'] + self.molecules['PFKb'] - totalPFK}")
            else: 
                #TODO: make assertion for "non-debug state"
                #assert 
                pass
                
            if self.verbose:
                print(f"Previous step values -- PFKb: {self.dictCountArrays['PFKb'][step-1]} PFKa: {self.dictCountArrays['PFKa'][step-1]} PFKi: {self.dictCountArrays['PFKi'][step-1]}")
                print(f"Current step values -- PFKb: {self.dictCountArrays['PFKb'][step]} PFKa: {self.dictCountArrays['PFKa'][step]} PFKi: {self.dictCountArrays['PFKi'][step]}")
                print(f"Difference -- PFKb: {self.dictCountArrays['PFKb'][step] - self.dictCountArrays['PFKb'][step-1]} PFKa: {self.dictCountArrays['PFKa'][step] - self.dictCountArrays['PFKa'][step-1]} PFKi: {self.dictCountArrays['PFKi'][step] - self.dictCountArrays['PFKi'][step-1]}")

        
        totalADPATP = self.dictCountArrays['ATP'][0] + self.dictCountArrays['ADP'][0] + self.dictCountArrays['PFKa'][0]
        if abs(self.molecules['ATP'] + self.molecules['ADP'] + self.molecules['PFKa'] - totalADPATP) > 1e-6:
            if self.debug:
                warnings.warn(f"On Step: {step} ATP-ADP mass Not Conserved! \n Correct Total: {totalADPATP}  Actual Total: {self.molecules['ATP'] + self.molecules['ADP'] + self.molecules['PFKa']}  Difference: {self.molecules['ATP'] + self.molecules['ADP'] + self.molecules['PFKa']- totalADPATP}")
            else: 
                #TODO: make assertion for "non-debug state"
                #assert 
                pass
            
            if self.verbose:
                # print(f"Previous step values -- ATP: {self.dictCountArrays['ATP'][step-1]} ADP: {self.dictCountArrays['ADP'][step-1]} ADP-PFKa: {self.dictCountArrays['PFKa'][step-1]}")
                # print(f"Current step values --  ATP: {self.dictCountArrays['ATP'][step]} ADP: {self.dictCountArrays['ADP'][step]} ADP-PFKa: {self.dictCountArrays['PFKa'][step]}")
                # print(f"Difference -- ATP: {self.dictCountArrays['ATP'][step] - self.dictCountArrays['ATP'][step-1]} ADP: {self.dictCountArrays['ADP'][step] - self.dictCountArrays['ADP'][step-1]} ADP-PFKa: {self.dictCountArrays['PFKa'][step]-self.dictCountArrays['PFKa'][step-1]}")
                print(f"Passive ATP Consumption: {self.passATPConsumption} Passive PEP Consumption: {self.passPEPConsumption} Passive PYR Consumption: {self.passPYRConsumption}")

        # Assert no concentration goes subatomic
        for key in self.dictCountArrays.keys():
            #assert (self.dictCountArrays[key][step] > approx1MoleculeInuM), f"On step {step}, {key} has gone into subatomic counts!!"
            if self.dictCountArrays[key][step] < ONE_COUNT_IN_MICROMOLAR and self.dictCountArrays[key][0] > ONE_COUNT_IN_MICROMOLAR:
                warnings.warn(f"On step {step}, {key} has gone into subatomic counts!!") 

        # Quasi-Steady State Validity Check:
        # E / (S + Km) << 1
        for key in self.dictQSSACheck.keys():
            # key is formatted: enzymeName_substrate thus: key.split("_")[0] extracts the enzyme
            enz, subs = key.split("_")[0], key.split("_")[1]

            qssa = self.molecules[enz] / (self.molecules[subs] + self.dictKM[f"{enz}_{subs}"])
            self.dictQSSACheck[key][step] = qssa > self.QSSAThreshold # Note: True indicates a violation of QSSA
            
            # Flag the QSSA violation to the threshold
            if self.dictQSSACheck[key][step] and not self.QSSAFlag[key]:
                warnings.warn(f"On step {step}, {enz} has violated the QSSA for {subs}!!") 
                self.QSSAFlag[key] = True


    def printStatements(self, step, rawRate):
        if step % 1000 == 0:
                self.percComplete = round((step*100 / self.time_totalSteps),1)
                printBlockMessage(f"Step: {step}. Simulation is {self.percComplete}% Complete.")
                print(f"Is Molecule Concentration 0? {self.zerosCheck}")
                print(f"Accumulated ATP-ADP error: {self.molecules['ATP'] + self.molecules['ADP'] + self.molecules['PFKa'] - (self.dictCountArrays['ATP'][0] + self.dictCountArrays['ADP'][0] + self.dictCountArrays['PFKa'][0])}")
                print(f"Molecule Turnover rates: { {key:self.dictRawStepTurnover[key][step] for key in self.dictRawStepTurnover.keys()} }")
                #print(rawRate)

    def run(self):
        printBlockMessage("Begin simulation...", ts =True)
        for step in range(1,self.time_totalSteps):       
            # Calculate raw concentration changes for current step.
            rawRate = self.rate()

            # Debugging: Save the raw turnover rate
            for key in rawRate.keys():
                self.dictRawStepTurnover[key][step] = rawRate[key]
                
            # Mass Conservation
            self.zerosCheck = { key:self.molecules[key] + rawRate[key] > ONE_COUNT_IN_MICROMOLAR for key in self.dictCountArrays.keys()}

            for key in self.dictCountArrays.keys():
                self.dictAdjustedStepTurnover[key][step] = rawRate[key] * self.zerosCheck[key]
                self.dictCountArrays[key][step] = self.molecules[key] + self.dictAdjustedStepTurnover[key][step]

            # Save step vmax:
            #self.saveStepVmax(step)

            # update current counts
            self.updateMolecules(step)
            
            # Check for purturbation:
            self.perturbation(step)

            if self.debug:
                self.checkAssertions(step)
            
            if self.verbose:
                self.printStatements(step, rawRate)

        printBlockMessage("Simulation Complete!", ts =True)
        return self.time, self.dictCountArrays, self.dictEnzymePercentMaxActivty, self.dictQSSACheck


if __name__ == '__main__':
    # visit http://127.0.0.1:8050/ in your web browser.
    from dash import Dash, html, dcc, Input, Output, State, dash_table, no_update
    from dash.dash_table.Format import Format, Scheme, Trim
    from dash.exceptions import PreventUpdate
    from vis import * #pplot, pplot_ioEnz, pplot_sat
    
    TABLE_STYLE_1 = {'width': '100px','maxWidth': '200px','minWidth': '100px',}

    app = Dash(__name__)

    app.config.suppress_callback_exceptions = True

    app.layout = html.Div([
        # The plot and time slider
        dcc.Graph(id='simulation-plot'),
        dcc.RangeSlider(id='x-step-slider',
            # min = min(),
            # max = max(time),
            # step=None,
            # value=[min(time), max(time)],
            # marks=None,
        ),
        html.Br(),
        html.Div([
            html.Div([ # Plot Options e
                html.P('Metabolite Concentration Plot: y-axis scale'),
                dcc.RadioItems(id = 'simulation-plot-metabolite-yaxis',
                    options = [
                        {'label':'Log-Scale' , 'value':'log'},
                        {'label':'Linear-Scale' , 'value':'linear'}
                        ],
                    value = 'log',
                    inline = True
                ),
                # Add Raw change from initial value 
                dcc.Checklist(id='delta-t-zero-checkbox',
                    options = [{'label':'Subtract t0 metabolite concentration from all values?', 'value':0}]
                ),
                dcc.Checklist(id='show-qsaa-checkbox',
                    options = [{'label':'Show where QSSA Violated?', 'value':0}]
                ),
            ], style = {'width': '40%', 'display':'inline-block'}),
            html.Div([ # Time Datatable
                dash_table.DataTable(id = 'sim-time-table', 
                    columns = ([
                        dict(id = 'simDurationSeconds', name = 'simDurationSeconds'),
                        dict(id = 'simStepsPerSecond', name = 'simStepsPerSecond')
                    ]),
                    data = [DICT_TIME],
                    editable = True,
                    style_data=TABLE_STYLE_1
                ),
            ], style = {'width': '60%', 'display':'inline-block'}),
            
            # # DICT_ENZ_SUBSTRATE
            # html.Div([ # Enz-Subs Pair Datatable
            #     dash_table.DataTable(id = 'ESP-table', 
            #         columns = ([
            #             dict(id = 'enz', name = 'enz'),
            #             dict(id = 'subs', name = 'subs'),
            #         ]),
            #         data = [dict(**{'Enzyme': key, 'Substrate': subs } ) for key in DICT_ENZ_SUBSTRATE.keys() for subs in DICT_ENZ_SUBSTRATE[key]],
            #         editable = False,
            #         style_data= TABLE_STYLE_1
            #     ),
            # ], style = {'width': '30%', 'display':'inline-block'}),
        ]),
        # K-cat value container/ k-cat button
        html.Button('Resimulate', id='resim-button', n_clicks=0),html.Button('Plot', id='plot-button', n_clicks=0),
        
        # Data Tables holding key Values for simulation
        html.Div([
            html.Div([ # Conc Datatable
                dash_table.DataTable(id = 'conc-val-table', 
                    columns = ([
                        dict(id = 'Name', name = 'Name'),
                        dict(id = 'Conc', name = 'Conc'),
                        dict(id = 'Type', name = 'Type',)
                    ]),
                    data = [dict(Name = key, **{'Conc': float( DICT_T0_CONC[key][0]), 'Type': DICT_T0_CONC[key][1] } ) for key in DICT_T0_CONC.keys()],
                    editable = True,
                    style_data= TABLE_STYLE_1
                ),
            ], style = {'width': '20%', 'display':'inline-block'}),
            html.Div([ # KM datatable
                dash_table.DataTable(id = 'km-val-table', 
                    columns = ([
                        dict(id = 'Name', name = 'Name'),
                        dict(id = 'KM', name = 'KM', format = Format(precision = 3, scheme = Scheme.exponent) )
                    ]),
                    data = [dict(Name = i, **{'KM': float( DICT_KM[i])}) for i in DICT_KM],
                    editable = True,
                    style_data=TABLE_STYLE_1
                )
            ], style = {'width': '20%', 'display':'inline-block'}),
            html.Div([ # k-val Datatable
                dash_table.DataTable(id = 'k-val-table', 
                    columns = ([
                        dict(id = 'Name', name = 'Name'),
                        dict(id = 'K', name = 'K', format = Format(precision = 3, scheme = Scheme.exponent) )
                    ]),
                    data = [dict(Name = i,**{'K': float(DICT_K[i])} ) for i in DICT_K],
                    editable = True,
                    style_data=TABLE_STYLE_1
                ),
            ], style = {'width': '20%', 'display':'inline-block'}),    
            html.Div([ # kcat-val Datatable
                dash_table.DataTable(id = 'kcat-val-table', 
                    columns = ([
                        dict(id = 'Name', name = 'Name'),
                        dict(id = 'Kcat', name = 'Kcat', format = Format(precision = 3, scheme = Scheme.exponent) )
                    ]),
                    data = [dict(Name = i,**{'Kcat': float(DICT_KCAT[i])} ) for i in DICT_KCAT],
                    editable = True,
                    style_data=TABLE_STYLE_1
                ),
            ], style = {'width': '20%', 'display':'inline-block'}),     
            html.Div([ # kreg-val Datatable
                dash_table.DataTable(id = 'kreg-val-table', 
                    columns = ([
                        dict(id = 'Name', name = 'Name'),
                        dict(id = 'Kreg', name = 'Kreg', format = Format(precision = 3, scheme = Scheme.exponent) )
                    ]),
                    data = [dict(Name = i,**{'Kreg': float(DICT_K_REGULATORY[i])} ) for i in DICT_K_REGULATORY],
                    editable = True,
                    style_data=TABLE_STYLE_1
                ),
            ], style = {'width': '20%', 'display':'inline-block'}),

        ]),
        html.Div([
            dcc.Dropdown(id = "ioEnz-enzName", placeholder="Select an Enzyme",
                options = [key for key in DICT_ENZ_SUBSTRATE.keys()],
            ),
            dcc.Dropdown(id = "ioEnz-x-axis", placeholder="Select a Substrate"),#, multi = True),
            dcc.Dropdown(id = "ioEnz-y-Axis", placeholder="Select Another Substrate or Modifier"),
            dcc.Dropdown(id = "ioEnz-y-Axis-type", placeholder="Select Modifer Type"),
            html.Button('Plot', id='plot-ioEnz-button', n_clicks=0),
            
        ]),
        html.Div([
            dcc.Graph(id='ioEnz-plot'),

        ]),
        
        # dcc.store 

        # store values 
        dcc.Store(id = 'simOut-time'),
        dcc.Store(id = 'simOut-stepSize'),
        dcc.Store(id = 'simOut-counts'),
        dcc.Store(id = 'simOut-rawStepTurnover'),
        dcc.Store(id = 'simOut-adjustedStepTurnover'),
        dcc.Store(id = 'simOut-enzymePercentSubstrateSaturation'),
        dcc.Store(id = 'simOut-enzymePercentMaxActivity'),
        dcc.Store(id = 'simOut-qssa'),

        # Inputs to Sim
        dcc.Store(id = 'sim-time', data = DICT_TIME),
        dcc.Store(id = 'protein-k-values'),
        dcc.Store(id = 'protein-km-values'),
        dcc.Store(id = 'protein-conc-values'),
        dcc.Store(id = 'protein-kcat-values'),
        dcc.Store(id = 'protein-kreg-values'),
    ])
    
    # Enzyme IO WIP
    @app.callback(
        Output('ioEnz-x-axis', 'options'),
        Input('ioEnz-enzName', 'value'),
        )
    def updateIOEnzXAxisOptions(data): 
        if data is None:
            raise PreventUpdate
        return [ o[0] for o in DICT_ENZ_SUBSTRATE[data]['subs']]
    # Enzyme IO WIP
    @app.callback(
        Output('ioEnz-y-Axis', 'options'),
        Input('ioEnz-enzName', 'value'),
        Input('ioEnz-x-axis', 'value'),
        )
    def updateIOEnzYAxisOptions(enzName, dataX): 
        if dataX is None or enzName is None:
            raise PreventUpdate
        return list(set([o[0] for key in DICT_ENZ_SUBSTRATE[enzName] for o in DICT_ENZ_SUBSTRATE[enzName][key] if o != dataX and key != 'prod']))

    @app.callback(
        Output('ioEnz-y-Axis-type', 'options'),
        Input('ioEnz-enzName', 'value'),
        Input('ioEnz-y-Axis', 'value'),
        )
    def updateIOEnzYAxisTypeOptions(enzName, dataY): 
        if enzName is None or dataY is None:
            raise PreventUpdate
        return list(set( o[1] for key in DICT_ENZ_SUBSTRATE[enzName] for o in DICT_ENZ_SUBSTRATE[enzName][key] if o[0] == dataY))


    @app.callback(
        Output('sim-time-table', 'data'),
        Input('sim-time-table', 'data_timestamp'),
        State('sim-time-table', 'data'))
    def updateTimeTable(timestamp, data): 
        if data is None:
            raise PreventUpdate
        return data
    @app.callback(
        Output('conc-val-table', 'data'),
        Input('conc-val-table', 'data_timestamp'),
        State('conc-val-table', 'data'))
    def updateConcTable(timestamp, data): 
        if data is None:
            raise PreventUpdate
        return data
    @app.callback(
        Output('km-val-table', 'data'),
        Input('km-val-table', 'data_timestamp'),
        State('km-val-table', 'data'))
    def updateKMTable(timestamp, data): 
        if data is None:
            raise PreventUpdate
        return data
    @app.callback(
        Output('k-val-table', 'data'),
        Input('k-val-table', 'data_timestamp'),
        State('k-val-table', 'data'))
    def updateKTable(timestamp, data): 
        if data is None:
            raise PreventUpdate
        return data
    @app.callback(
        Output('kcat-val-table', 'data'),
        Input('kcat-val-table', 'data_timestamp'),
        State('kcat-val-table', 'data'))
    def updateKcatTable(timestamp, data): 
        if data is None:
            raise PreventUpdate
        return data
    @app.callback(
        Output('kreg-val-table', 'data'),
        Input('kreg-val-table', 'data_timestamp'),
        State('kreg-val-table', 'data'))
    def updateKregTable(timestamp, data): 
        if data is None:
            raise PreventUpdate
        return data
    @app.callback(
        Output('protein-conc-values', 'data'),
        Input('conc-val-table', 'data'),)
    def updateConc(data):
        if data is None:
            raise PreventUpdate
        return data
    @app.callback(
        Output('protein-km-values', 'data'),
        Input('km-val-table', 'data'))
    def updateKM(data):
        if data is None:
            raise PreventUpdate
        return data
    @app.callback(
        Output('protein-k-values', 'data'),
        Input('k-val-table', 'data'))
    def updateK(data):
        if data is None:
            raise PreventUpdate
        return data
    @app.callback(
        Output('protein-kcat-values', 'data'),
        Input('kcat-val-table', 'data'))
    def updateKcat(data):
        if data is None:
            raise PreventUpdate
        return data  
    @app.callback(
        Output('sim-time', 'data'),
        Input('sim-time-table', 'data'))
    def updateSimTime(data):
        if data is None:
            raise PreventUpdate
        return data   
    @app.callback(
        Output('protein-kreg-values', 'data'),
        Input('kreg-val-table', 'data'))
    def updateKcat(data):
        if data is None:
            raise PreventUpdate
        return data


    @app.callback(
        Output('simOut-time', 'data'),
#        Output('simOut-stepSize', 'data'),
        Output('simOut-counts', 'data'),
        # Output('simOut-rawStepTurnover', 'data'),
        # Output('simOut-adjustedStepTurnover', 'data'),
        # Output('simOut-enzymePercentSubstrateSaturation', 'data'),
        Output('simOut-enzymePercentMaxActivity', 'data'),
        Output('simOut-qssa', 'data'),
        
        Input('resim-button', 'n_clicks'),
        State('sim-time', 'data'),
        State('protein-conc-values', 'data'),
        State('protein-kcat-values', 'data'),
        State('protein-k-values', 'data'),
        State('protein-km-values', 'data'),
        State('protein-kreg-values', 'data'), 
        prevent_initial_call = True
    )
    def updateSim(resimbutton, time, conc, kcat, k, km,kreg):
        # Convert data to correct format:
        timeInput = {key:int(time[0][key]) for key in time[0].keys()}
        concInput = {conc[i]['Name']:(float(conc[i]['Conc']), conc[i]['Type']) for i in range(len(conc))}
        kcatInput = {kcat[i]['Name']: float(kcat[i]['Kcat']) for i in range(len(kcat))}
        kInput =  {k[i]['Name']: float(k[i]['K']) for i in range(len(k))}
        kmInput =  {km[i]['Name']: float(km[i]['KM']) for i in range(len(km))}
        kregInput =  {kreg[i]['Name']: float(kreg[i]['Kreg']) for i in range(len(kreg))}

        sim = Glycolysis(timeInput, concInput, kcatInput, kInput, kmInput, kregInput)
        outT, outC, outPercVmax, qssa = sim.run()
        return outT, outC, outPercVmax, qssa


    @app.callback(
        Output('simulation-plot', 'figure'),
        Input('plot-button', 'n_clicks'),
        Input('simulation-plot-metabolite-yaxis', 'value'),
        Input('delta-t-zero-checkbox', 'value'),
        Input('show-qsaa-checkbox', 'value'),
        #Input('x-step-slider', 'value'),
        State('simOut-time', 'data'),
        #State('simOut-stepSize', 'data'),
        State('simOut-counts', 'data'),
        State('simOut-qssa', 'data'),
        # State('simOut-rawStepTurnover', 'data'),
        # State('simOut-adjustedStepTurnover', 'data'),
        # State('simOut-enzymePercentSubstrateSaturation', 'data'),
        # State('simOut-enzymePercentMaxActivity', 'data'),
        prevent_initial_call = True
        )
    def update_figure(plot_button, metaboliteYAxis, deltaTZero, showQSSA, time, count, qssa):
        # convert data to the Plotly_Dynamics format
        fig=pplot(
                Time = time,
                metaboliteCountData = count,
                qssa = qssa,
                yAxisScaleMetabolite=metaboliteYAxis,
                deltaT0Scale=deltaTZero,
                showQSSA = showQSSA,
                excludedMolecules=EXCLUDE_FROM_PLOT,
                )
        return fig

    @app.callback( # ioEnz plot
        Output('ioEnz-plot', 'figure'),
        Input('plot-ioEnz-button', 'n_clicks'),
        
        # Get enzyme, X and Y axis
        State('ioEnz-enzName', 'value'),
        State('ioEnz-x-axis', 'value'), # X is always sat for now
        State('ioEnz-y-Axis', 'value'),
        State('ioEnz-y-Axis-type', 'value'),
        
        # Sim stuff
        State('simOut-counts', 'data'),
        State('simOut-time', 'data'),
        State('protein-km-values', 'data'),
        State('protein-kreg-values', 'data'), 
        prevent_initial_call = True
        )
    def update_ioEnz_figure(plot_button,enzName, XName, YName, YType , count, time, km, kreg):
        # convert data to the Plotly_Dynamics format
       
        # Convert data to correct format:
        kmInput =  {km[i]['Name']: float(km[i]['KM']) for i in range(len(km))}
        kregInput =  {kreg[i]['Name']: float(kreg[i]['Kreg']) for i in range(len(kreg))}
        
        fig=pplot_sat_sim(enzName, XName, YName, YType, time, count,
                            #kcatInput, kInput,
                             kmInput, kregInput)
        return fig


    app.run_server(debug=True)







