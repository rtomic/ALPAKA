HELP := $(word 1, $(TMP_MODULES))
TMP_MODULES := $(wordlist 2,  $(words $(TMP_MODULES)), $(TMP_MODULES))

############################################################################################################
####Module Objects
TMP_OBJ = \
 ErrorFactory \
 LooError \
 RandomPointsError \
 CError \
 #L2Error  RandomPointsError#Object-Files   Error  LooError
####SubModules
TMP_MOD =  #Falls Unterordner mit anderen Modules vorhanden, die da rein, sonst leer
############################################################################################################


OBJECTS += $(patsubst %, $(HELP)/%.o, $(TMP_OBJ))

SUB_MODS := $(patsubst %, $(HELP)/%, $(TMP_MOD))

MODULES += $(SUB_MODS)
TMP_MODULES := $(SUB_MODS) $(TMP_MODULES) 

include $(patsubst %, %/module.mk, $(SUB_MODS))
