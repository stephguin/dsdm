# Makefile generique Multi-plateforme pour construire une appli avec :

#   - tous les .c des repertoires  */C   ou C
#
# ------------------------------------------------------------------------------------------
# Variables generales 
# ------------------------------------------------------------------------------------------
HOST    = $(shell hostname)
CUR_DIR = $(shell pwd) 
SYSTEM  = $(shell uname -s)
RDIST_DIR = delami
#
# ------------------------------------------------------------------------------------------
#   - Quelques options de compilateur et editeur de lien pour memoire
#
#   - Quelques options interressantes de g77  ( faire info -f g77 )
#       -Wunused
#       -Wuninitialized
#       -Wall ( les 2 precedentes )
#       -fno_backslash  ( \ n'est pas un caractere special )
#       -fdollar-ok     ( accepte le $ dans les noms )
#       -fugly-assumed  ( x(1) devient x(*) )
#
#       -finit_local_zero
#       -fno_underscoring
#
#
#   - Quelques options de ld
#       -d ou -dc ou -dp pour compatibilite avec d'autres linker
#          permet l'allocation d'espace aux common et evite des
#          warning etranges ( ???? )
#
#
# ------------------------------------------------------------------------------------------
# Parametres de compilation et d'edition de liens suivant le system
# ------------------------------------------------------------------------------------------
#
#   ---------------------------------------- Linux
#
ifeq  "$(SYSTEM)"  "Linux" 
 CC  = gcc
 F77 = gfortran
 MAKE = make
#
 CFLAGS  = -O
##FFLAGS  = -Wimplicit -fugly -fno-underscoring -fno-backslash -fdollar-ok -Iincl
 FFLAGS  = -fno-underscoring -fdollar-ok -Iincl -fd-lines-as-comments
# FFLAGS  = -Wall -fno-underscoring -fno-backslash -fdollar-ok -Iincl
 LDFLAGS = -d
#
 CFLAGSD  = -g 
# FFLAGSD  = -g -D__TRACE__ -Wimplicit -fno-underscoring -fno-backslash -fdollar-ok -Iincl
 FFLAGSD  = -g -fno-underscoring -fdollar-ok -Iincl -fd-lines-as-code
 LDFLAGSD = -g -d
endif
#
#   ---------------------------------------- HP-UX
#
ifeq "$(SYSTEM)"  "HP-UX"
 CC  = gcc
 F77 = f77
 MAKE = gmake
#
 CFLAGS  = -O
 FFLAGS  = +O3 -Iincl
 LDFLAGS =
#
 CFLAGSD  = -g
 FFLAGSD  = -g
 LDFLAGSD = -g
endif
#
#
#   ---------------------------------------- OSF1
#
#   Pour digital, l'option de compil -pad_source supprime les messages
#   d'info zarbi style "Hollerith constant"
#
ifeq "$(SYSTEM)"  "OSF1"
 CC  = gcc
 F77 = f77
 MAKE = gmake
#
 CFLAGS  = -O
#FFLAGS  = +O3 -Iincl -pad_source
 FFLAGS  = +O3 -Iincl
 LDFLAGS =
#
 CFLAGSD  = -g
 FFLAGSD  = -g
 LDFLAGSD = -g
endif
#
# -------------------------------------------------------------------------------------------
# Listes des sources des objets et des includes
# -------------------------------------------------------------------------------------------
#
# Liste des sources ( wildcard GNU)
SRCF  = $(wildcard */ftn/*.f)
SRCC  = $(wildcard */C/*.c)
#
# Liste des objets ( syntaxe GNU)
OBJF  = $(SRCF:.f=.o)
OBJC  = $(SRCC:.c=.o)
OBJFD  = $(SRCF:.f=.oo)
OBJCD  = $(SRCC:.c=.oo)
#
# Liste des fichiers include ou header
INCL = $(wildcard incl/*.h)
INHH = $(INCL:.h=.hh)
#
# -------------------------------------------------------------------------------------------
#  Editeur de lien suivant les types de sources 
#  et les librairies associees aux langages.
#--------------------------------------------------------------------------------------------
#
#    - Pour memoire, quelques infos sur les librairies
# 
#        La commande suivante donne le nom de la lib gcc chargee automatiquement
#        gcc -print-libgcc-file-name
#        /usr/lib/gcc-lib/i486-linux/egcs-2.91.60/libgcc.a
#
#        pour f77 on doit probablement linker
#          /usr/lib/gcc-lib/i486-linux/egcs-2.91.60/libg2c.a
#          /usr/lib/gcc-lib/i486-linux/egcs-2.91.60/libgcc.a
#        en effet faire 
#          g77 -print-file-name=libg2c.a
#          g77 -print-libgcc-file-name
#    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    - ? - ? - ? - ? - ? - ? - ? - ? - ? - ? - ? - ? - ? - ? - ? - ? - ? - ? - ? - ? -
#    -  
#    - Problemes restant a regler
#       - des warning aussi abondants que bizarres sur les commons !!!! 
#         -d=-dc=-dp       une option de ld            ne marche pas
#         une version plus recente du compilo devrait resoudre le PB
#         (voir compil de castem par Panpan a Versailles)
#         pour l'instant on fait des grep -v
#       - d'une maniere generale, l'utilisation generique de LD !!
#    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#
LIBLG = 

ifeq "$(SRCF)" ""
  LD = $(CC)
endif
ifeq "$(SRCC)" ""
#
# ---------------------------------------------------------------------------------------
# avec g77 et toutes ces options , on n'arrive pas a eliminer les messages des commons
#
#  LD = $(F77)
#  LDFLAGS := -Wl,-g, -dc   pas de message de common mais pb no such file or directory
#                           doit planter avant les commons !!!
#  LDFLAGS := -Wl,-g,-dc    MESSAGES sur les common et no such
#  LDFLAGS := -Wl,-g  -dc   idem
#  LDFLAGS := -dc           idem
#
#  LDFLAGS :=  -Wl,g,-dc,    pas de mess common mais pb nosuch
#  LDFLAGS :=  -Wl,          idem !
#
#
# ---------------- avec ce ld il manque   _start, il n'y a plus de message sur les commons
#    mais ils seiraient surement apres
#    LD = ld
#    LIBLG = -L/usr/lib/gcc-lib/i486-linux/egcs-2.91.60  -lgcc -lg2c -lstdc++
#    LDFLAGS = -d
# -----------------------------------------------------------------------------------------
#
endif
ifneq "$(SRCC)" ""
 ifneq "$(SRCF)" ""
   LD = ld
   ifeq  "$(SYSTEM)"  "Linux" 
     LIBLG = -L/usr/lib/gcc-lib/i486-linux/egcs-2.91.60  -lgcc -lg2c -lstdc++
   endif
 endif
endif

# -------------------------------------------------------------------------------------------
# Eventuellement des Librairies
#                                                  !!! a customizer avec soin !!!
# -------------------------------------------------------------------------------------------
#
LIB = /home/sg/dsdm/interpret/lib_new.1 \
      /home/sg/dsdm/gesdyn/lib_util.1 \
      /home/sg/dsdm/divers/lib_div
#
LID = /home/sg/dsdm/interpret/lib_new_d.1 \
      /home/sg/dsdm/gesdyn/lib_util_d.1 \
      /home/sg/dsdm/divers/lib_div_d
#
# Dans un premier temps on ne met pas lib_ester, tout le code a ete
# rapatrie en local
#LIB = /u/lmtutils/gesdyn/lib_util    \
#      /u/lmtutils/interpret/lib_new.1 
#      /u/ester/pp2/lib_ester
#
# -------------------------------------------------------------------------------------------
# Nom de l'executable ou de la librairie
#                                                  !!! a customizer avec soin !!!
# -------------------------------------------------------------------------------------------
#
EXE      = dsdm
EXED     = dsdm_d
LIB_EXE  =
LIB_EXED = 
#
# -------------------------------------------------------------------------------------------
default : n 
#
.PHONY : all dezero clean info n d h help corton vougeot romanee bowmore
#
# -------------------------------------------------------------------------------------------
# differents racourcis pour les cibles
#
all : $(EXE) $(EXED) $(LIB_EXE) $(LIB_EXED)
n   : $(EXE)
d   : $(EXED)
lib : $(LIB_EXE)
h   : help
#
# -------------------------------------------------------------------------------------------
# - Le help
# -------------------------------------------------------------------------------------------
help :
	@echo '------------------------------- help --------------------------------------------'
	@echo '---------------- Makefile generique pour presque tout faire ---------------------'
	@echo '------------------------------- help --------------------------------------------'
	@echo ""
	@echo "    make, make h ou make help  permet d'obtenir cette aide"
	@echo "    make n           pour obtenir un executable en compilation normale"
	@echo "    make d           . . . . . . . . . . . . .  en mode debug"
	@echo "    make all         pour faire les 2 versions, normale et debug"
	@echo "    make lib         pour faire une librairie plutot qu'un executable"
	@echo "    make info        pour connaitre les options, les sources ..."
	@echo "    make clean       pour supprimer les relogeables et les executables"
	@echo "    make dezero      pour repartir de zero en recompilant tout"
	@echo ""
	@echo "    make corton ou vougeot ou romanee ou bowmore  "
	@echo "                     pour mettre a jour la version sur la machine nommee"
#
# -------------------------------------------------------------------------------------------
# Regle de contruction de l'executable NORMAL 
# -------------------------------------------------------------------------------------------
#
$(EXE) :  $(INHH) $(OBJC) $(OBJF)
	  @echo "On va construire l'executable $(EXE) "
	  $(F77) -o $@ $(OBJF) $(LIB) 2>&1
#
# -------------------------------------------------------------------------------------------
# Regle de contruction de la librairie NORMALE 
# -------------------------------------------------------------------------------------------
#
$(LIB_EXE) :  $(INHH) $(OBJC) $(OBJF)
	  @echo "On va construire la librairie $(LIB_EXE) "
	  ar -cr $@   $(OBJC) $(OBJF)
#
# -------------------------------------------------------------------------------------------
# Regle de contruction de l'executable DEBUG  
# -------------------------------------------------------------------------------------------
# 
# edition de lien avec g77 en DUR !!! et supression des warnings !!!
#
$(EXED) :  $(INHH) $(OBJCD) $(OBJFD)
	@echo "On va construire l'executable DEBUG $(EXED) "
	$(F77) -g -o $@ $(OBJFD) $(LID) 2>&1 
	
#	| grep -v 'multiple common' \
#	                           | grep -v 'previous common is here'     \
#	                           | grep -v 'overridden by definition' \
#	                           | grep -v 'overridden by larger common' \
#	                           | grep -v 'larger common is here' \
#	                           | grep -v 'defined here' \
#	                           | grep -v 'overriding common' \
#	                           | grep -v 'overriding smaller common' \
#	                           | grep -v 'common is here' 
#	  g77 -g -o $@ $(OBJFD) std_gest/C/*.o $(LIB) 2>&1   | grep -v 'multiple common'
#  
# regle generale a faire marcher des que possible !!!
#
#	  $(LD) -o $@ $(LDFLAGSD)  $(OBJCD) $(OBJFD) $(LIB) $(LIBLG)
#
# -------------------------------------------------------------------------------------------
# Regle de contruction de la librairie DEBUG  
# -------------------------------------------------------------------------------------------
#
$(LIB_EXED) :  $(INHH) $(OBJCD) $(OBJFD)
	  @echo "On va construire la librairie $(LIB_EXED) "
	  ar -cr $@   $(OBJCD) $(OBJFD)
#
# -------------------------------------------------------------------------------------------
# Certains fichiers include ont ils ete modifies, si oui on cherche les sources qui y font
# reference et on les touche.
# -------------------------------------------------------------------------------------------
#               - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#               - Quelques notations pour memoire
#
#                 $(<F)  $(<D) la partie fichier et repertoire de la premiere
#                              dependance que l'on note $<
#               - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#
%.hh : %.h
	@echo 'le header $@ a ete modifie'
ifneq "$(SRCF)" ""
	-@grep -le ^[^Cc].*$(<F) */ftn/*.f  | xargs touch
endif
ifneq "$(SRCC)" ""
	-@grep -le ^#.*$(<F) C/*.c  | xargs touch
endif
	@touch $@        
 
# -------------------------------------------------------------------------------------------
#               - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#               - Quelques notations pour memoire
#
#                 $@  la cible                         --- Variables automatiques - p 92
#                 $+  la liste de toutes les dependance
#                 $<  la pemiere dependance
#                 $?  liste des dependances qui on changees
#                 $*  le stem (radical), partie de la cible correspondant a une regle implicite
#                        si %.o et toto.o -> toto
#               - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# -------------------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------------------
# Regles implicites de compilation fortran et C
# -------------------------------------------------------------------------------------------
#
# Si fichier *.F, g77 appelle automatiquement cpp
#
%.o : %.f
#	@echo 'compile d un .f'
	$(F77) $(FFLAGS) -c $*.f   -o $@

%.o : %.c
#	@echo 'compile d un .c'
	$(CC)  $(CFLAGS) -c $<   -o $@
#
%.oo : %.f
#	@echo 'compile debug d un .f'
	$(F77) $(FFLAGSD) -c $*.f   -o $@

%.oo : %.c
#	@echo 'compile debug d un .c'
	$(CC)  $(CFLAGSD) -c $<   -o $@
#
#
# -----------------------------------------------------------------------------------------------
# Label pour une mise a jour entre systemes, appel a rdist
#    ATTENTION il faut mettre a jour le Rdist_file_gene
# -----------------------------------------------------------------------------------------------
#
vougeot corton romanee bowmore :
	@echo "Copy de $(HOST) vers $@"
	rdist -f Rdist_file_gene  $(RDIST_DIR)_$(HOST)_$@
#
#
# -----------------------------------------------------------------------------------------------
# Label pour repartir de zero, et tout recompiler
# -----------------------------------------------------------------------------------------------
#
dezero : 
	@touch $(SRCF) $(SRCC)
	@$(MAKE) all
#
# -----------------------------------------------------------------------------------------------
# Label pour obtenir des infos sur les parametres et variables
# -------------------------------------------------------------------------------------------
#
info :
	@echo "------------------------------------------ Infos generales "
	@echo $(CUR_DIR)
	@echo $(SYSTEM)
	@echo $(HOST)
	@echo "----------------------------------------- Fichiers include "
	@echo $(INCL)
	@echo "------------------------------------------------- Sources C"
	@echo $(SRCC)
	@echo "------------------------------------------- Sources FORTRAN"
	@echo $(SRCF)
	@echo "-------------------------------------------------- Objets C"
	@echo $(OBJC)
	@echo "-------------------------------------------- Objets FORTRAN"
	@echo $(OBJF)
	@echo "--------------------------------------- Infos compilation C"
	@echo "Compilateur : " $(CC)
	@echo "Options std : " $(CFLAGS)
	@echo "Options dbg : " $(CFLAGSD)
	@echo "--------------------------------- Infos compilation FORTRAN"
	@echo "Compilateur : " $(F77)
	@echo "Options std : " $(FFLAGS)
	@echo "Options dbg : " $(FFLAGSD)
	@echo "------------------------------------------ Editeur de liens"
	@echo "Editeur de lien : " $(LD)
#
# -------------------------------------------------------------------------------------------
# Label pour nettoyer les objets et les executables
# -------------------------------------------------------------------------------------------
#
clean :
	find . -name '*.o' -exec rm {} \;
	find . -name '*.oo' -exec rm {} \;
	rm -f $(EXE) $(EXED) $(LIB_EXE) $(LIB_EXED)


















