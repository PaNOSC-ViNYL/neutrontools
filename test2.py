from SimEx.Parameters.IonMatterInteractorParameters import IonMatterInteractorParameters
from SimEx.Calculators.TNSAIonMatterInteractor import TNSAIonMatterInteractor

myparams = IonMatterInteractorParameters(ion_name='proton', neutron_weight=1.e4)
#List of parameters: energy_bin, neutron_weight, ibeam_radius, target_length, target_density,
                    #source_dump, xsec_file
                    #ion_name
myparams.xsec_file = 'D_D_-_3He_n.txt'

mysource = TNSAIonMatterInteractor(parameters=myparams, input_path='Data/0010.sdf', output_path='Data/NeutronData.h5')

mysource.run()

mysource.saveH5()

