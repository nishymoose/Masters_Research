import machupX as MX

import numpy as np
from math import exp, pow, sqrt, cos, sin, tan, atan, asin, acos, atan2
import math
import matplotlib.pyplot as plt
import json
import numpy.matlib
from scipy import linalg
import scipy.linalg as spla

# imported to replicate figure 1.8.19
import math  as m
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import chi2, loglaplace
from labellines import labelLine, labelLines
from random import randint
from operator import itemgetter
# imported to replicate figure 1.8.19


class Masters_Research:
    def __init__(self, test_condition):
        
        
        
        if test_condition == "IGE":
            # initializing dictionary in ground effect
            self.data = {
                "CG" : [0,0,0],
                "weight" : 100.0,
                "reference" : {
                    "area" : 225.0,
                    "longitudinal_length" : 6.5,
                    "lateral_length" : 36.74
                },
                "airfoils" : {
                    "NACA_0010" : {
                        "type" : "linear",
                        "aL0" : 0.0,
                        "CLa" : 6.28318530717959,
                        "CmL0" : 0.0,
                        "Cma" : 0.00,
                        "CD0" : 0.00,
                        "CD1" : 0.0,
                        "CD2" : 0.0,
                        "CL_max" : 1.4
                    }
                },
                "wings" : {
                    "main_wing" : {
                        "ID" : 1,
                        "side" : "both",
                        "is_main" : True,
                        "connect_to" : {
                            "location" : "root",
                            "dz" : -1.837
                        },
                        "semispan" : 18.37,
                        "airfoil" : "NACA_0010",
                        "twist" : [[0.0, 0.0], [1.0, 0.0]],
                        "chord" : [[0.0, 8.75], [1.0, 3.5]],
                        "grid" : {
                            "N" : 40
                        }
                    },
                    "ground_effect_wing" : {
                        "ID" : 2,
                        "side" : "both",
                        "is_main" : False,
                        "connect_to" : {
                            "location" : "root",
                            "dz" : 1.837
                        },
                        "semispan" : 18.37,
                        "airfoil" : "NACA_0010",
                        "twist" : [[0.0, 0.0], [1.0, 0.0]],
                        "chord" : [[0.0, 8.75], [1.0, 3.5]],
                        "grid" : {
                            "N" : 40
                        }
                    }
                }
            }
            
            # Writes Dictionairy for Wing IN GROUND EFFECT(IGE) to json to be used in MachupX calculations
            with open("traditional_airplane_case_one.json", "w") as write_file:
                json.dump(self.data, write_file, indent = 2)
                
                
        elif test_condition == "OGE":
            # initializing dictionary out of ground effect
            self.data = {
                "CG" : [0,0,0],
                "weight" : 100.0,
                "reference" : {
                    "area" : 225.0,
                    "longitudinal_length" : 6.5,
                    "lateral_length" : 36.74
                },
                "airfoils" : {
                    "NACA_0010" : {
                        "type" : "linear",
                        "aL0" : 0.0,
                        "CLa" : 6.28318530717959,
                        "CmL0" : 0.0,
                        "Cma" : 0.00,
                        "CD0" : 0.00,
                        "CD1" : 0.0,
                        "CD2" : 0.0,
                        "CL_max" : 1.4
                    }
                },
                "wings" : {
                    "main_wing" : {
                        "ID" : 1,
                        "side" : "both",
                        "is_main" : True,
                        "semispan" : 18.37,
                        "airfoil" : "NACA_0010",
                        "twist" : [[0.0, 0.0], [1.0, 0.0]],
                        "chord" : [[0.0, 8.75], [1.0, 3.5]],
                        "grid" : {
                            "N" : 40
                        }
                    }
                }
            }
            
            # Writes Dictionairy for Wing NOT IN GROUND EFFECT(OGE) to json to be used in MachupX calculations
            with open("traditional_airplane_case_one.json", "w") as write_file:
                json.dump(self.data, write_file, indent = 2)
            
        
        elif test_condition == "elliptic":
            # initializing dictionary out of ground effect
            self.data = {
                "CG" : [0,0,0],
                "weight" : 100.0,
                "reference" : {
                    "area" : 225.0,
                    "longitudinal_length" : 6.5,
                    "lateral_length" : 36.74
                },
                "airfoils" : {
                    "NACA_0010" : {
                        "type" : "linear",
                        "aL0" : 0.0,
                        "CLa" : 6.28318530717959,
                        "CmL0" : 0.0,
                        "Cma" : 0.00,
                        "CD0" : 0.00,
                        "CD1" : 0.0,
                        "CD2" : 0.0,
                        "CL_max" : 1.4
                    }
                },
                "wings" : {
                    "main_wing" : {
                        "ID" : 1,
                        "side" : "both",
                        "is_main" : True,
                        "semispan" : 18.37,
                        "airfoil" : "NACA_0010",
                        "twist" : [[0.0, 2.0], [1.0, 2.0]],
                        "chord" : ["elliptic", 8.0],
                        "grid" : {
                            "N" : 40
                        }
                    }
                }
            }
            
            # Writes Dictionairy for Wing NOT IN GROUND EFFECT(OGE) to json to be used in MachupX calculations
            with open("traditional_airplane_case_one.json", "w") as write_file:
                json.dump(self.data, write_file, indent = 2)
            
        else:
            print("Incorrect Input")
            
        
        self.test_case = test_condition
        
        # Specify input file and Create global myScene
        self.input_file = "traditional_input_case_one.json"
        self.my_scene = MX.Scene(self.input_file)
        
        ################# FORCES #################
        self.my_forces = self.my_scene.solve_forces(dimensional=False, non_dimensional=True, report_by_segment=True, verbose=False)
        
        self.CL_total_main = self.my_forces["traditional_airplane"]["inviscid"]["CL"]["main_wing_left"] + self.my_forces["traditional_airplane"]["inviscid"]["CL"]["main_wing_right"]
        self.CD_total_main= self.my_forces["traditional_airplane"]["inviscid"]["CD"]["main_wing_left"] + self.my_forces["traditional_airplane"]["inviscid"]["CD"]["main_wing_right"]

        ################# DISTRIBUTIONS #################
        self.distributions = self.my_scene.distributions(verbose=False)
        
        #print("CD_i distributions", self.distributions["traditional_airplane"]["main_wing_left"]["CD_i"])
        #print()
        
        # CONSTANTS
        self.chord_length = self.data["wings"]["main_wing"]["chord"]
        self.root_chord_list = self.chord_length[0]
        self.tip_chord_list = self.chord_length[1]
        
        self.c_root = self.root_chord_list[1]
        if test_condition != "elliptic":
            self.c_tip = self.tip_chord_list[1]
        self.semispan = self.data["wings"]["main_wing"] ["semispan"]
        self.b = self.semispan*2
        self.area = self.data["reference"]["area"]
        self.RA = (self.b*self.b)/self.area
        
        
        
        if test_condition == "IGE":
            self.main_wing_height_above_ground = self.data["wings"]["main_wing"]["connect_to"]["dz"]
            self.ground_effect_wing_height_above_ground = self.data["wings"]["ground_effect_wing"]["connect_to"]["dz"]
            
        
        if test_condition != "elliptic":
            self.RT = self.c_tip/self.c_root

        self.CLa = self.data["airfoils"]["NACA_0010"]["CLa"]
        
        #TODO: Do I want to prompt the user to input his value?
        self.CL = 0.5
        
    # Displays wing plot
    def Display_scene(self):
        self.my_scene.display_wireframe(show_legend=True)
        
    # Displays wing plot
    def Display_planform(self):
        self.my_scene.display_planform()
        
    # Returns the angle of attack required to hit specified CL for given planform  (DOES NOT APPLY THEM TO CALCULATIONS, ONLY RETURNS A VALUE)
    def Get_AoA_At_Target_CL(self, target_CL):
        return self.my_scene.target_CL(CL = target_CL)     
        
    
    
    
    # Applies angle of attack needed to hit target CL to plaform for calculations
    def Apply_AoA_At_Target_CL(self, target_CL):
        
        
        Applied_AoA = self.Get_AoA_At_Target_CL(target_CL)
        print()
        print("Applied_AoA: ", Applied_AoA)
        print()
        
        def Extract(lst):
            return [item[0] for item in lst]
        
        def set_AoA_Main_Wing(lst, x):
            return [[item[0], item[1]+x]for item in lst]
        
        def set_AoA_Ground_Effect_Wing(lst, x):
            return [[item[0], item[1]+x]for item in lst]
        
        base_twist_list_Main_Wing = self.data["wings"]["main_wing"] ["twist"]
        if self.test_case == "IGE":
            base_twist_list_Ground_Effect_Wing = self.data["wings"]["ground_effect_wing"] ["twist"]
            
        
        print()
        print("CL Total: ",json.dumps(self.my_forces["traditional_airplane"]["inviscid"]["CL"]["total"], indent=4))
        print("CD Total: ",json.dumps(self.my_forces["traditional_airplane"]["inviscid"]["CD"]["total"], indent=4))
        print()
        
        with open("traditional_airplane_case_one.json", "w") as write_file:
            if self.test_case == "IGE":
                self.data["wings"]["main_wing"]["twist"]= set_AoA_Main_Wing(base_twist_list_Main_Wing, Applied_AoA)
                self.data["wings"]["ground_effect_wing"]["twist"]= set_AoA_Ground_Effect_Wing(base_twist_list_Ground_Effect_Wing, -Applied_AoA)
            
            elif self.test_case == "OGE":
                self.data["wings"]["main_wing"]["twist"] = set_AoA_Main_Wing(base_twist_list_Main_Wing, Applied_AoA)
                #self.my_scene.display_wireframe(show_legend=True)
            
            elif self.test_case == "elliptic":
                self.my_scene.display_wireframe(show_legend=True)
                
            else:
                print("Incorrect Input")
            json.dump(self.data, write_file, indent = 2)
        
        self.my_scene = MX.Scene(self.input_file)
        #self.my_scene.display_wireframe(show_legend=True)
        #print("twist after AoA sweep", self.data["wings"]["main_wing"]["twist"])
        self.my_forces = self.my_scene.solve_forces(dimensional=False, non_dimensional=True, report_by_segment=True, verbose=False)
        
        # Prints overall Lift and Drag, CL should be zero if done correctly
        print()
        print("CL Total: ",json.dumps(self.my_forces["traditional_airplane"]["inviscid"]["CL"]["total"], indent=4))
        print("CD Total: ",json.dumps(self.my_forces["traditional_airplane"]["inviscid"]["CD"]["total"], indent=4))
        print()
        
        
        self.distributions = self.my_scene.distributions(verbose=False)
        #self.Plot_CL_Distribution()
        
        
    # Reversing a list using reversed()
    def Reverse(self, lst):
        return [ele for ele in reversed(lst)]
    
    
    # Need the span fraction for each node to complete twist optimization
    def get_span_frac(self):
        self.span_frac = self.distributions["traditional_airplane"]["main_wing_left"]["span_frac"]
        return self.span_frac
    
    
    def get_span_frac_left(self):
        self.span_frac_left= self.distributions["traditional_airplane"]["main_wing_left"]["span_frac"]
        return self.span_frac_left
    
    def get_span_frac_right(self):
        self.span_frac_right = self.distributions["traditional_airplane"]["main_wing_right"]["span_frac"]
        return self.span_frac_right
    
    # Eq. 1.8.33 in Phillips for chord length for a tapered wing
    def calculate_c_z(self, z_spanwise, RT):
        return (2*self.b/(self.RA*(1+RT)))*(1 - (1 - RT)*(np.absolute((2*z_spanwise)/self.b)))


    # Eq. 1.8.33 in Phillips for chord length for a tapered wing
    def calculate_c_theta(self, theta, RT):
        return (2*self.b/(self.RA*(1+RT)))*(1 - (1 - RT)*(np.absolute(m.cos(theta))))
    
    # change of variables for the spansise coordinate in liftine-line theory, pg 52 in Phillips
    # Eq. 1.8.3 in Phillips
    def calculate_theta(self, z_spanwise):
        return m.acos(-2*z_spanwise/self.b)
    
    
    # washout distribution function
    # Eq. 1.8.42 in Phillips ()
    def calculate_omega(self, theta, c_theta):
        return 1 - (m.sin(theta)/(c_theta/self.c_root))


    # Optimum total washout to minimize induced drag
    # Eq. 1.8.42 in Phillips
    def calculate_Cap_omega(self):
        return 4*self.b/(m.pi*self.RA*self.CLa*self.c_root)
    
    
    # Optimum total washout to minimize induced drag
    # Eq. 1.8.43 in Phillips
    def linear_taper_calculate_Cap_omega(self, RT):
        return (2*(1 + RT)*self.CL)/(m.pi*self.CLa)
    
    
    # washout distribution function
    # Eq. 1.8.43 in Phillips ()
    def linear_taper_calculate_omega(self, theta, RT):
        return 1 - (m.sin(theta)/(1 - (1 - RT)*(np.absolute(m.cos(theta)))))
    
    
    # washout distribution function
    # Eq. 1.8.42 in Phillips ()
    def calculate_omega(self, theta, c_theta):
        return 1 - (m.sin(theta)/(c_theta/self.c_root))
    
    
    # Relicates figure 1.8.19 in Phillips textbook
    def Display_Optimum_Washout_Distribution_Plot(self):    
        
# =============================================================================
#         # Eq. 1.8.33 in Phillips for chord length for a tapered wing
#         def calculate_c_z(self, z_spanwise):
#             return (2*self.b/(self.RA*(1+RT)))*(1 - (1 - RT)*(np.absolute((2*z_spanwise)/self.b)))
# 
#         # Eq. 1.8.33 in Phillips for chord length for a tapered wing
#         def calculate_c_theta(self, theta):
#             return (2*self.b/(self.RA*(1+RT)))*(1 - (1 - RT)*(np.absolute(m.cos(theta))))
# 
#         # change of variables for the spansise coordinate in liftine-line theory, pg 52 in Phillips
#         # Eq. 1.8.3 in Phillips
#         def calculate_theta(self, z_spanwise):
#             return m.acos(-2*z_spanwise/self.b)
# 
#         # washout distribution function
#         # Eq. 1.8.42 in Phillips ()
#         def calculate_omega(self, theta, c_theta):
#             return 1 - (m.sin(theta)/(c_theta/self.c_root))
# 
#         # Optimum total washout to minimize induced drag
#         # Eq. 1.8.42 in Phillips
#         def calculate_Cap_omega(self):
#             return 4*self.b/(m.pi*self.RA*self.CLa*self.c_root)
# 
#         # washout distribution function
#         # Eq. 1.8.43 in Phillips ()
#         def linear_taper_calculate_omega(self, theta):
#             return 1 - (m.sin(theta)/(1 - (1 - RT)*(np.absolute(m.cos(theta)))))
# 
#         # Optimum total washout to minimize induced drag
#         # Eq. 1.8.43 in Phillips
#         def linear_taper_calculate_Cap_omega(self):
#             return (2*(1 + RT)*self.CL)/(m.pi*self.CLa)
# =============================================================================
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #print(self.get_span_frac())
        for i in np.arange(0,1.1,0.1):
            RT = i

            # START ============================================================================= STEP 1:
            # section span_frac will give me the location of each node along the span in a value from 0 to 1
            # multiply this value by the semispan to get the spanwise location z
            span_frac = self.get_span_frac()
           
            z_spanwise = [element * self.semispan for element in span_frac]
            # END ============================================================================= STEP 1
            
            # START ============================================================================= STEP 2: 
            # change of variables for spanwise coordinate
            change_of_variables = []
            list_length = len(span_frac)
            
            for i in range(list_length):
               new_theta = self.calculate_theta(z_spanwise[i])
               change_of_variables.append(new_theta)
               
            # theta at each spanwise location is now equal to change_of_variables, it is an array of length 40
            theta = change_of_variables
            # END ============================================================================= STEP 2: 
            
            # START ============================================================================= STEP 3:  
            # are already given the chord length
            c_z_array = []
            
            c_theta_array = []
            
            for i in range(list_length):
               c_z_array.append(self.calculate_c_z(z_spanwise[i], RT))
               c_theta_array.append(self.calculate_c_theta(theta[i], RT) )
            # END ============================================================================= STEP 3:
            
            # START ============================================================================= STEP 4:  
            omega_opt_array = []
            for i in range(list_length):
               omega_opt_array.append(self.linear_taper_calculate_omega(theta[i], RT))
            # END ============================================================================= STEP 4:
                
            # START ============================================================================= STEP 4b:  
            total_opt_twist = self.linear_taper_calculate_Cap_omega(RT)
            # this is given in radians, convert to degrees
            total_opt_twist_degrees = m.degrees(total_opt_twist)
            # END ============================================================================= STEP 4b:
                
            # START ============================================================================= STEP 5:  
            # this value is the value that will be inserted into the json file to twist the wing
            applied_twist = []
            for i in range(list_length):
               apply_twist_at_node = omega_opt_array[i] * total_opt_twist_degrees
               applied_twist.append(apply_twist_at_node)
            # END ============================================================================= STEP 5:
               
                
            # START ============================================================================= STEP 6:  
            # Replicate plot in Figure 1.8.19 in Phillips
            z_over_b = []
            for i in range(list_length):
                z_over_b_at_node = z_spanwise[i]/self.b
                z_over_b.append(z_over_b_at_node)
                
            x = z_over_b
            y = omega_opt_array
            
            # All lines are the same color
            #ax.plot(x,y, 'k', lw=2, label = RT)
            
            # All lines are different colors
            RT_rounded = RT.round(decimals=2)
        
            ax.plot(x,y, lw=2, label = RT_rounded)
            
            ax.hlines(y=0.0, xmin=0.0, xmax=1.0, color='k', linestyle = 'dashed')
            
            ax.vlines(x=0.1, ymin=0.0, ymax=0.25, color='k')
            
            ax.text(0.025, 0.25, 'elliptic planform (\u03C9 = 0)')
            
            ax.set_ylim(-1.0,1.0)
            ax.set_xlim(0.0, 0.5)
            
        # Show plot outside of for loop to show all lines on same graph       
        ax.minorticks_on()
        ax.tick_params(which='major', length=10, width=2, direction='inout')
        ax.tick_params(which='minor', length=5, width=2, direction='in')
        ax.grid(which='both')

        labelLines(plt.gca().get_lines(), align=False, zorder=2.5)

        plt.title("Figure 1.8.19")

        #txt="I need the caption to be present a little below X-axis"

        ax.set_xlabel('z/b \n Optimum washout distribution for wings of linear taper \n in production of minimum induced drag, as defined in Eq. 1.8.43 in Phillips')
        plt.ylabel("\u03C9", rotation=0)
        
        plt.show()
        
    

    def replicate_figure_3_from_hunsaker_paper(self):
        test = 0
        
        
        
        
    # Calculates the optimal twist distribution for the given planform and applies this twist to the self.data dictionairy
    def Apply_Optimal_Twist_Distribution(self):
        
        span_frac = self.get_span_frac()
        list_length = len(span_frac)
        
        z_spanwise = [element * self.semispan for element in span_frac]
        #print("z_spanwise: ", z_spanwise)
        #print()
        # Optimum total washout to minimize induced drag
        # Eq. 1.8.43 in Phillips
        total_opt_twist = self.linear_taper_calculate_Cap_omega(self.RT)
        # this is given in radians, convert to degrees
        total_opt_twist_degrees = m.degrees(total_opt_twist)
        #print("total_opt_twist_degrees: ", total_opt_twist_degrees)
        #print()
        
        
        # change of variables for spanwise coordinate
        change_of_variables = []
        
        for i in range(list_length):
           new_theta = self.calculate_theta(z_spanwise[i])
           change_of_variables.append(new_theta)
           
        # theta at each spanwise location is now equal to change_of_variables, it is an array of length 40
        theta = change_of_variables
        #print()
        #print("theta: ", theta)
        #print()
        
        omega_opt_array = []
        for i in range(list_length):
           omega_opt_array.append(self.linear_taper_calculate_omega(theta[i], self.RT))
        
        
        #print()
        #print("omega_opt_array: ", omega_opt_array)
        #print()
        # calculates the twist value to be applies at each span fraction location
        applied_twist = []
        for i in range(list_length):
           apply_twist_at_node = omega_opt_array[i] * total_opt_twist_degrees
           applied_twist.append(apply_twist_at_node)
        
        #print()
        #print("applied_twist: ", applied_twist)
        #print()
        #print("applied_twist: ", applied_twist)
        #print()
        # If i want to create a lists of lists to apply the twist to the wing planform in the dictionary
        # Concerned that there are actually 42 specified locations instead of 40 because I added specifications
        # for values at 0.0 and 1.0
        
        
        l = []
        l.append([1.0, total_opt_twist_degrees])  
        for i in range(list_length):
           l.append([span_frac[i], applied_twist[i]])
           
        l.append([0.0, 0.0])
        #print(l)
        twist_list = self.Reverse(l)
        
        
        
        
        applied_twist_ground_effect_wing = [element * -1 for element in applied_twist]
        
        
        
        l_ground_effect_wing = []
        l_ground_effect_wing.append([1.0, -total_opt_twist_degrees])  
        for i in range(list_length):
           l_ground_effect_wing.append([span_frac[i], applied_twist_ground_effect_wing[i]])
           
        l_ground_effect_wing.append([0.0, 0.0])
        #print(l)
        twist_list_ground_effect_wing = self.Reverse(l_ground_effect_wing)
        
        
        # concerned that twist_list has 42 elements instead of 40 and that that will mess with my values
        
        
        if self.test_case == "IGE":
            # initializing dictionary in ground effect
            self.data = {
                "CG" : [0,0,0],
                "weight" : 100.0,
                "reference" : {
                    "area" : 225.0,
                    "longitudinal_length" : 6.5,
                    "lateral_length" : 36.74
                },
                "airfoils" : {
                    "NACA_0010" : {
                        "type" : "linear",
                        "aL0" : 0.0,
                        "CLa" : 6.28318530717959,
                        "CmL0" : 0.0,
                        "Cma" : 0.00,
                        "CD0" : 0.00,
                        "CD1" : 0.0,
                        "CD2" : 0.0,
                        "CL_max" : 1.4
                    }
                },
                "wings" : {
                    "main_wing" : {
                        "ID" : 1,
                        "side" : "both",
                        "is_main" : True,
                        "connect_to" : {
                            "location" : "root",
                            "dz" : -1.837
                        },
                        "semispan" : 18.37,
                        "airfoil" : "NACA_0010",
                        "twist" : twist_list,
                        "chord" : [[0.0, 8.75], [1.0, 3.5]],
                        "grid" : {
                            "N" : 40
                        }
                    },
                    "ground_effect_wing" : {
                        "ID" : 2,
                        "side" : "both",
                        "is_main" : False,
                        "connect_to" : {
                            "location" : "root",
                            "dz" : 1.837
                        },
                        "semispan" : 18.37,
                        "airfoil" : "NACA_0010",
                        "twist" : twist_list_ground_effect_wing,
                        "chord" : [[0.0, 8.75], [1.0, 3.5]],
                        "grid" : {
                            "N" : 40
                        }
                    }
                }
            }
            
            # Writes Dictionairy for Wing IN GROUND EFFECT(IGE) to json to be used in MachupX calculations
            with open("traditional_airplane_case_one.json", "w") as write_file:
                json.dump(self.data, write_file, indent = 2)
                
                
        elif self.test_case == "OGE":
            # initializing dictionary out of ground effect
            self.data = {
                "CG" : [0,0,0],
                "weight" : 100.0,
                "reference" : {
                    "area" : 225.0,
                    "longitudinal_length" : 6.5,
                    "lateral_length" : 36.74
                },
                "airfoils" : {
                    "NACA_0010" : {
                        "type" : "linear",
                        "aL0" : 0.0,
                        "CLa" : 6.28318530717959,
                        "CmL0" : 0.0,
                        "Cma" : 0.00,
                        "CD0" : 0.00,
                        "CD1" : 0.0,
                        "CD2" : 0.0,
                        "CL_max" : 1.4
                    }
                },
                "wings" : {
                    "main_wing" : {
                        "ID" : 1,
                        "side" : "both",
                        "is_main" : True,
                        "semispan" : 18.37,
                        "airfoil" : "NACA_0010",
                        "twist" : twist_list,
                        "chord" : [[0.0, 8.75], [1.0, 3.5]],
                        "grid" : {
                            "N" : 40
                        }
                    }
                }
            }
            
            # Writes Dictionairy for Wing NOT IN GROUND EFFECT(OGE) to json to be used in MachupX calculations
            with open("traditional_airplane_case_one.json", "w") as write_file:
                json.dump(self.data, write_file, indent = 2)
            
        else:
            print("Incorrect Input")
        
        #print("twist", self.data["wings"]["main_wing"] ["twist"])
        # Specify input file and Create global myScene
        self.input_file = "traditional_input_case_one.json"
        self.my_scene = MX.Scene(self.input_file)
        
        ################# FORCES #################
        self.my_forces = self.my_scene.solve_forces(dimensional=False, non_dimensional=True, report_by_segment=True, verbose=False)
        
        with open('twist_distribution_results_with_dict.json', 'w', encoding='utf-8') as f:
            json.dump(self.my_forces, f, ensure_ascii=False, indent=4)
            
        self.CL_total_main = self.my_forces["traditional_airplane"]["inviscid"]["CL"]["main_wing_left"] + self.my_forces["traditional_airplane"]["inviscid"]["CL"]["main_wing_right"]
        self.CD_total_main= self.my_forces["traditional_airplane"]["inviscid"]["CD"]["main_wing_left"] + self.my_forces["traditional_airplane"]["inviscid"]["CD"]["main_wing_right"]
       
        
        #print("Total CL: ", self.CL_total_main)
        #print("Total CD: ", self.CD_total_main)
        
        # IF WING IS IN GROUND EFFECT, C_TOTAL SHOULD ALWAYS BE EQUAL TO 0
        print("-------------------------------------------")
        print()
        print("CL Total: ",json.dumps(self.my_forces["traditional_airplane"]["inviscid"]["CL"]["total"], indent=4))
        print("CD Total: ",json.dumps(self.my_forces["traditional_airplane"]["inviscid"]["CD"]["total"], indent=4))
        print()
        print("-------------------------------------------")
        
        print("-------------------------------------------")
        print()
        print("CL Main Wing: ",json.dumps(self.my_forces["traditional_airplane"]["inviscid"]["CL"]["main_wing_left"] + self.my_forces["traditional_airplane"]["inviscid"]["CL"]["main_wing_right"], indent=4))
        print("CD Main Wing: ",json.dumps(self.my_forces["traditional_airplane"]["inviscid"]["CD"]["main_wing_left"] + self.my_forces["traditional_airplane"]["inviscid"]["CD"]["main_wing_right"], indent=4))
        print()
        print("-------------------------------------------")
        
        ################# DISTRIBUTIONS #################
        self.distributions = self.my_scene.distributions(verbose=False)   
        #print("twist", self.distributions["traditional_airplane"]["main_wing_left"] ["twist"])
        with open('twist_distribution_results_with_dict.json', 'w', encoding='utf-8') as f:
            json.dump(self.distributions, f, ensure_ascii=False, indent=4)
    
    # Plots the section CL values for each node of the given planform
    def Plot_CL_Distribution(self):
        ################# DISTRIBUTIONS #################
        left_section_CL = self.distributions["traditional_airplane"]["main_wing_left"] ["section_CL"]
        right_section_CL = self.distributions["traditional_airplane"]["main_wing_right"] ["section_CL"]
        
        combined_section_CL = left_section_CL + right_section_CL

        left_span_frac = self.get_span_frac_left()
        
        
        right_span_frac = self.get_span_frac_right()
        
        multiplied_left_span_frac = [element * -1 for element in left_span_frac]
        
        total_span_frac = multiplied_left_span_frac + right_span_frac
        
        fig,ax = plt.subplots(1)
        
        x = total_span_frac
        y = combined_section_CL

        # plot the data
        ax.plot(x,y)
    
    # performs a sweep over a range of angles of attack specified by the user and returns the Lift and Drag coefficients
    def AoA_sweep(self, range_low, range_high, step):
        
        def Extract(lst):
            return [item[0] for item in lst]
        
        def adjust_AoA_Main_wing(lst, x):
            return [[item[0], item[1]+x]for item in lst]
        
        def adjust_AoA_Ground_Effect_Wing(lst, x):
            return [[item[0], item[1]+x]for item in lst]
        
        base_twist_list_Main_Wing = self.data["wings"]["main_wing"] ["twist"]
        if self.test_case == "IGE":
            base_twist_list_Ground_Effect_Wing = self.data["wings"]["ground_effect_wing"] ["twist"]
        
        #AoA_Swept_twist_list = adjust_AoA(base_twist_list, 1)
        #print("twist after AoA sweep", AoA_Swept_twist_list)
        for x in range(range_low, range_high, step):
            with open("traditional_airplane_case_one.json", "w") as write_file:
                if self.test_case == "IGE":
                    self.data["wings"]["main_wing"]["twist"]= adjust_AoA_Main_wing(base_twist_list_Main_Wing, x)
                    self.data["wings"]["ground_effect_wing"]["twist"]= adjust_AoA_Ground_Effect_Wing(base_twist_list_Ground_Effect_Wing, -x)
                    
                elif self.test_case == "OGE":
                    self.data["wings"]["main_wing"]["twist"] = adjust_AoA_Main_wing(base_twist_list_Main_Wing, x)
                    self.my_scene.display_wireframe(show_legend=True)
                
                elif self.test_case == "elliptic":
                    self.my_scene.display_wireframe(show_legend=True)
                    
                else:
                    print("Incorrect Input")
                json.dump(self.data, write_file, indent = 2)
            self.my_scene = MX.Scene(self.input_file)
            self.my_scene.display_wireframe(show_legend=True)
            #print("twist after AoA sweep", self.data["wings"]["main_wing"]["twist"])
            self.my_forces = self.my_scene.solve_forces(dimensional=False, non_dimensional=True, report_by_segment=True, verbose=False)
            
            # Prints overall Lift and Drag, CL should be zero if done correctly
            #print()
            
            # IF WING IS IN GROUND EFFECT, C_TOTAL SHOULD ALWAYS BE EQUAL TO 0
            print("-------------------------------------------")
            print()
            print("CL Total: ",json.dumps(self.my_forces["traditional_airplane"]["inviscid"]["CL"]["total"], indent=4))
            print("CD Total: ",json.dumps(self.my_forces["traditional_airplane"]["inviscid"]["CD"]["total"], indent=4))
            print()
            print("-------------------------------------------")
            
            print("-------------------------------------------")
            print()
            print("CL Main Wing: ",json.dumps(self.my_forces["traditional_airplane"]["inviscid"]["CL"]["main_wing_left"] + self.my_forces["traditional_airplane"]["inviscid"]["CL"]["main_wing_right"], indent=4))
            print("CD Main Wing: ",json.dumps(self.my_forces["traditional_airplane"]["inviscid"]["CD"]["main_wing_left"] + self.my_forces["traditional_airplane"]["inviscid"]["CD"]["main_wing_right"], indent=4))
            print()
            print("-------------------------------------------")
            
            print("AoA: ", x)
            
            
            self.distributions = self.my_scene.distributions(verbose=False)
            #self.Plot_CL_Distribution()
            
    
        
    # If in ground effect, sets the height of the wing above the ground
    def set_h_over_b(self,h_over_b_in):
        if self.test_case == "IGE":
            
            # These are the current values for height above ground before applying the new height
            print()
            print("-----------------------------------------------")
            print("-----------------------------------------------")
            print("Previous main wing height above ground: ", np.absolute(self.main_wing_height_above_ground))
            print("Previous ground effect wing height below ground: ", self.ground_effect_wing_height_above_ground)
            h_over_b_ratio = np.absolute(self.main_wing_height_above_ground/self.b)
            print("Previous (h/b) ratio of height above ground to wingspan: ", h_over_b_ratio)
            print()
            #print("-----------------------------------------------")
            print("------------SETTING NEW WING HEIGHT------------")
            #print("-----------------------------------------------")
            print()
            
            
            
            
            with open("traditional_airplane_case_one.json", "w") as write_file:
                if self.test_case == "IGE":
                    self.data["wings"]["main_wing"]["connect_to"]["dz"] = -1*(h_over_b_in*self.b)
                    self.data["wings"]["ground_effect_wing"]["connect_to"]["dz"] = -self.main_wing_height_above_ground
                    
                    self.main_wing_height_above_ground = self.data["wings"]["main_wing"]["connect_to"]["dz"]
                    self.ground_effect_wing_height_above_ground = self.data["wings"]["ground_effect_wing"]["connect_to"]["dz"]
                
                else:
                    print("Incorrect Input")
                    
                json.dump(self.data, write_file, indent = 2)
            self.my_scene = MX.Scene(self.input_file)
            self.my_scene.display_wireframe(show_legend=True)
            
            # need to set the location of dz in data to the value for it to show up I want it to
            
            # These are the new values for height above ground after applying the new height
            
            
            h_over_b_ratio_out = np.absolute(self.main_wing_height_above_ground/self.b)
            print("New main wing height above ground: ", np.absolute(self.main_wing_height_above_ground))
            print("New ground effect wing height below ground: ", self.ground_effect_wing_height_above_ground)
            print("New (h/b) ratio of height above ground to wingspan: ", h_over_b_ratio_out)
            print("-----------------------------------------------")
            print("-----------------------------------------------")
            print()
            
        else:
            print("Invalid entry, height above ground only callable in ground effect")
        
    # Returns the height above the ground of the main wing and the h/b ratio
    def get_height_above_ground_and_h_over_b(self):
 
        if self.test_case == "IGE":
            print("Main wing height above ground: ", np.absolute(self.main_wing_height_above_ground))
            print("Ground effect wing height below ground: ", self.ground_effect_wing_height_above_ground)
            h_over_b_ratio = np.absolute(self.main_wing_height_above_ground/self.b)
            print("(h/b) ratio of height above ground to wingspan: ", h_over_b_ratio)
        else:
            print("Invalid entry, height above ground only callable in ground effect")
            
    
    # Performs a sweep over range of values specified for the user for the wing above the ground in ground effect
    def h_over_b_sweep(self, range_low, range_high, step_size):
        
        for x in np.arange(range_low, range_high, step_size):
            self.set_h_over_b(x)
            
            self.my_scene = MX.Scene(self.input_file)
            self.my_scene.display_wireframe(show_legend=True)
            
    
        
    def display_heat_map(self):
        test = 0
    
    # Calculates the optimal twist distribution for the given planform (DOES NOT APPLY THESE TO SELF.DATA FOR CALCULATIONS)
    def get_applied_twist_values(self):
        span_frac = self.get_span_frac()
        list_length = len(span_frac)
        
        
        z_spanwise_reversed = [element * self.semispan for element in span_frac]
        # currently, the list is given to me with z_spanwise_reversed[0] being equal to the location at the tip
        # this reverses the order so that z_spanwise_reversed[0] is equal to the first node after the root
        
        
        # Reversing a list using reversed()
        def Reverse(lst):
            return [ele for ele in reversed(lst)]
        
        z_spanwise = Reverse(z_spanwise_reversed)
        
        
        # Optimum total washout to minimize induced drag
        # Eq. 1.8.43 in Phillips
        total_opt_twist = self.linear_taper_calculate_Cap_omega(self.RT)
        # this is given in radians, convert to degrees
        total_opt_twist_degrees = m.degrees(total_opt_twist)
        
        
        # change of variables for spanwise coordinate
        change_of_variables = []
        
        for i in range(list_length):
           new_theta = self.calculate_theta(z_spanwise[i])
           change_of_variables.append(new_theta)
           
        # theta at each spanwise location is now equal to change_of_variables, it is an array of length 40
        theta = change_of_variables
        
        
        omega_opt_array = []
        for i in range(list_length):
           omega_opt_array.append(self.linear_taper_calculate_omega(theta[i], self.RT))
        
        
        
        # calculates the twist value to be applies at each span fraction location
        applied_twist = []
        for i in range(list_length):
           apply_twist_at_node = omega_opt_array[i] * total_opt_twist_degrees
           applied_twist.append(apply_twist_at_node)
           
        # If i want to create a lists of lists to apply the twist to the wing planform in the dictionary
        # Concerned that there are actually 42 specified locations instead of 40 because I added specifications
        # for values at 0.0 and 1.0
        span_frac_correct_direction = self.Reverse(span_frac)  
        
        l = []
        l.append([0.0, 0.0])
        for i in range(list_length):
           l.append([span_frac_correct_direction[i], applied_twist[i]])
           
        l.append([1.0, total_opt_twist_degrees])  
        #print(l)
        twist_list = l
        # concerned that twist_list has 42 elements instead of 40 and that that will mess with my values
           
        return applied_twist, twist_list
    
    
    # Generatesa 2D plot with lines for each chord length, and colors them based on the amount of twist
    # given the optimal twist distribution for the planform
    def Generate_Optimal_Twist_Heat_Map(self):
        span_frac = self.get_span_frac()
        list_length = len(span_frac)
        
        
        z_spanwise_reversed = [element * self.semispan for element in span_frac]
        z_spanwise_left_wing = [ -x for x in z_spanwise_reversed]
        
        z_spanwise_left_wing.insert(0, -self.semispan)
        z_spanwise_left_wing.append(0.0)
        # currently, the list is given to me with z_spanwise_reversed[0] being equal to the location at the tip
        # this reverses the order so that z_spanwise_reversed[0] is equal to the first node after the root
        
        
        # Reversing a list using reversed()
        def Reverse(lst):
            return [ele for ele in reversed(lst)]
        
        z_spanwise_right_wing = Reverse(z_spanwise_reversed)

        
        z_spanwise_right_wing.insert(0, 0.0)
        z_spanwise_right_wing.append(self.semispan)
        #print("z_spanwise_left_wing", z_spanwise_left_wing)
        ##print("z_spanwise_right_wing", z_spanwise_right_wing)
        
        full_main_wing_x = z_spanwise_left_wing
        
        full_main_wing_x.extend(z_spanwise_right_wing)
        #print("full_main_wing_x", full_main_wing_x)
        #print("full_main_wing_x length", len(full_main_wing_x))

        y_magnitude = []
        list_length = len(z_spanwise_right_wing)
        
        for i in range(list_length-2):
           new_y_magnitude = self.calculate_c_z(z_spanwise_right_wing[i], self.RT)
           y_magnitude.append(new_y_magnitude)
          
        y_mag_left_wing = Reverse(y_magnitude)
        
        
        # for right wing
        y_mag_right_wing = y_magnitude
        #y_mag_right_wing = [ -x for x in y_mag_right_wing]
        y_mag_right_wing.insert(0, self.c_root)
        y_mag_right_wing.append(self.c_tip)
        
        
        #print("y_mag_left_wing", y_mag_left_wing)
        #print("y_mag_left_wing", len(y_mag_left_wing))
        y_mag_left_wing.insert(0, self.c_tip)
        y_mag_left_wing.append(self.c_root)
        #print("y_mag_left_wing", y_mag_left_wing)
        #print("y_mag_left_wing", len(y_mag_left_wing))
        #print("y_mag_right_wing", y_mag_right_wing)
        
        y_mag_left_wing.extend(y_mag_right_wing)
        full_main_wing_y_mag = y_mag_left_wing
        #print("full_main_wing_y_mag", full_main_wing_y_mag)
        #print("full_main_wing_y_mag length", len(full_main_wing_y_mag))
        
        x_coord = full_main_wing_x
        y_coord_end_point = [element * 0.25 for element in full_main_wing_y_mag]
        y_coord_start_point_positive = [element * 0.75 for element in full_main_wing_y_mag]
        y_coord_start_point = [ -x for x in y_coord_start_point_positive]
        
        
        #print("x_coord: ", x_coord)
        #print("y_coord_start_point: ", y_coord_start_point)
        #print("y_coord_end_point: ", y_coord_end_point)
        
        
        #x_coord_end_point = [a + b for a, b in product(listB, listA)]
        
        
        data1 = np.array(x_coord).flatten()
        data2 = np.array(y_coord_start_point).flatten()
        data3 = np.array(y_coord_end_point).flatten()
        
        x = data1
        y = np.concatenate([data2[:,None],data3[:,None]], axis=1)
        #print(y[0,0])
        
        # This plots just the chord lengths of the wing in black and not a heat map
# =============================================================================
#         #print('x=', x,'y=',y)
#         plt.plot((x,x),([i for (i,j) in y], [j for (i,j) in y]),c='black')
#         #plt.plot(x, [i for (i,j) in y], 'rs', markersize = 4)
#         #plt.plot(x, [j for (i,j) in y], 'bo', markersize = 4)
#         x1,x2,y1,y2 = plt.axis()  
#         plt.axis((x1,x2,-15,15))
#         plt.show()
# =============================================================================
        
        applied_twist, twist_list_side_one = self.get_applied_twist_values()
        vmin = min(applied_twist)
        vmax = max(applied_twist)
        
        list_length = len(data1)
        
        # Reversing a list using reversed()
        def Reverse(lst):
            return [ele for ele in reversed(lst)]
        
        twist_list_side_two = Reverse(twist_list_side_one)
        
        twist_list = twist_list_side_two + twist_list_side_one
        segments_z =[]
        #print(len(twist_list))
        na = np.array(twist_list)
        #print(na[:,1])
        both_sides_applied_twist = na[:,1]
        #print(len(both_sides_applied_twist))
        for i in range(list_length-2):
            segments_z.append(((x[i],y[i,0]),(x[i],y[i,1]),both_sides_applied_twist[i]))
        
        #print(segments_z)
        #[((4, 20), (1, 4), 16), ((16, 6), (10, 19), 12), ((19, 8), (1, 5), 3), ((14, 20), (3, 12), 10), ((17, 20), (11, 6), 11), ((18, 14), (17, 5), 0), ((6, 17), (7, 15), 16), ((9, 0), (5, 3), 10), ((2, 12), (11, 17), 13), ((1, 20), (6, 6), 18), ((6, 12), (2, 11), 1), ((18, 0), (18, 5), 13)]

        
        
        
        #segments_z = [((x,y(0)),(x,y(1)),applied_twist) for _ in range(data1)]
        
        #[j for (i,j) in y]
        norm = plt.Normalize(vmin, vmax)
        cm = plt.cm.rainbow
        sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
        for p1, p2, z in segments_z:
            x, y = zip(p1,p2)
            
            plt.plot(x, y, color=cm(norm(z)))
        
        # draw the colorbar, note that we pass explicitly the ScalarMappable
        plt.colorbar(sm)
        
        # I'm done, I'll show the results,
        # you probably want to add labels to the axes and the colorbar.
        x1,x2,y1,y2 = plt.axis()  
        plt.axis((x1,x2,-15,15))
        plt.show()
        
        
        
    
    
        
            







p1 = Masters_Research("IGE")

#VERIFIED CORRECT
#p1.Display_scene()

#VERIFIED CORRECT
#p1.Display_planform()

#VERIFIED CORRECT
#p1.Display_Optimum_Washout_Distribution_Plot()

#VERIFIED CORRECT
#p1.Apply_Optimal_Twist_Distribution()


#VERIFIED CORRECT
#p1.Generate_Optimal_Twist_Heat_Map()


#p1.get_height_above_ground_and_h_over_b()



#p1.set_h_over_b(1)


p1.h_over_b_sweep(0,1,0.1)

#p1.Plot_CL_Distribution()




#print(p1.Get_AoA_At_Target_CL(0.5))


#p1.Apply_AoA_At_Target_CL(0.5)
#p1.Display_scene()

# Input in format for rnage of sweep with size step
# (range_low, range_high, step)
#p1.AoA_sweep(-1,6,1)

#p1.Plot_CL_Distribution()

#print(p1.Display_Optimum_Washout_Distribution_Plot())
#print("twist", p1.data["wings"]["main_wing"] ["twist"])

