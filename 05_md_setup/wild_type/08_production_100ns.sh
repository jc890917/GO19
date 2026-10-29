#!/bin/bash
#SBATCH --job-name=prod_md
#SBATCH --partition=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=4-00:00:00
#SBATCH --output=prod_%j.log
#SBATCH --error=prod_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jccarrillo@uach.mx
# ============================================================================
#  Step 08 — Production MD, 100 ns, réplica independiente, GPU
#
#  Uso:
#     sbatch 08_production_100ns.sh <REP>
#  o mediante job array:
#     sbatch --array=1-3 08_production_100ns.sh
#
#  <REP> = número de réplica (1, 2, 3, ...). Si no se pasa como argumento,
#  se toma de $SLURM_ARRAY_TASK_ID (para uso con --array).
#  El --job-name y --output de arriba se pueden sobreescribir por réplica
#  pasándolos en la línea de sbatch (ver lanzar_3_replicas.sh) — las
#  opciones de línea de comandos de sbatch tienen prioridad sobre las
#  #SBATCH del script.
#
#  Cada réplica vive en su propia carpeta prod_rep<REP>/ y arranca desde
#  equil/npt4.gro con velocidades regeneradas (semilla distinta por réplica)
#  para que las 3 trayectorias sean estadísticamente independientes, no
#  copias del mismo checkpoint.
#
#  --time=4-00:00:00 alcanzaba para 50 ns; para 100 ns puede que necesites
#  más pared de tiempo o dejar que el job se corte y se vuelva a encolar
#  (el script detecta prod.cpt y continúa automáticamente — ver más abajo).
#  Pregunta en Yoltla si hay límite máximo de --time en la partición gpus.
# ============================================================================

source ./config.sh
load_gromacs

# ── Determinar número de réplica ────────────────────────────────────────────
REP="${1:-$SLURM_ARRAY_TASK_ID}"
if [[ -z "$REP" ]]; then
    echo "[ERROR] No se especificó número de réplica."
    echo "        Uso: sbatch 08_production_100ns.sh <REP>"
    echo "        o:   sbatch --array=1-3 08_production_100ns.sh"
    exit 1
fi

# ── Parámetros de esta corrida (100 ns, override de config.sh) ─────────────
PRODUCTION_NS=100
DT=0.002                                             # ps, 2 fs
PRODUCTION_STEPS=$(awk -v ns="$PRODUCTION_NS" -v dt="$DT" 'BEGIN{printf "%d", (ns*1000)/dt}')

PRODDIR="prod_rep${REP}"
mkdir -p "$PRODDIR"

echo ""
echo "══════════════════════════════════════════════════"
echo "  Step 08 — Production MD (${PRODUCTION_NS} ns) — Réplica ${REP}"
echo "  Job ID : $SLURM_JOB_ID"
echo "  Node   : $SLURMD_NODENAME"
echo "  Dir    : $PRODDIR"
echo "══════════════════════════════════════════════════"

# ── Check input ───────────────────────────────────────────────────────────
[[ ! -f "equil/npt4.gro" ]] && {
    echo "[ERROR] equil/npt4.gro not found — run Step 07 first."; exit 1; }

NCPU="${SLURM_CPUS_PER_TASK:-10}"

# ── Semilla reproducible pero distinta por réplica ─────────────────────────
GEN_SEED=$(( 100003 + REP * 7919 ))

# ── Write production MDP (solo si no existe, para resumes estables) ────────
if [[ ! -f "${PRODDIR}/prod.mdp" ]]; then
    cat > "${PRODDIR}/prod.mdp" << EOF
; Production MD — ${PRODUCTION_NS} ns, NPT, ${TEMP} K, no restraints
; Réplica ${REP} — velocidades regeneradas con semilla propia (gen_seed=${GEN_SEED})
integrator          = md
nsteps              = ${PRODUCTION_STEPS}   ; ${PRODUCTION_NS} ns
dt                  = ${DT}                 ; 2 fs timestep
; Output frequency
nstxout-compressed  = 5000     ; trajectory every 5000 steps (10 ps)
compressed-x-grps   = Solute   ; solute only: ~45 MB per 100 ns instead of ~5 GB
nstenergy           = 5000     ; energy log every 10 ps
nstlog              = 5000     ; md log every 10 ps
nstxout             = 0        ; no uncompressed positions
nstvout             = 0        ; no velocities
; Bonds
continuation        = no                    ; velocidades nuevas, no del cpt de equil
constraint_algorithm= lincs
constraints         = h-bonds
lincs_iter          = 1
lincs_order         = 4
; Neighbor searching
cutoff-scheme       = Verlet
nstlist             = 10
ns_type             = grid
; Electrostatics (PME)
coulombtype         = PME
rcoulomb            = 1.0
pme_order           = 4
fourierspacing      = 0.16
; Van der Waals
rvdw                = 1.0
; Temperature — V-rescale thermostat
tcoupl              = V-rescale
tc-grps             = System
tau_t               = 0.1
ref_t               = ${TEMP}
; Pressure — Parrinello-Rahman (correct NPT ensemble for production)
pcoupl              = Parrinello-Rahman
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = ${PRESSURE}
compressibility     = 4.5e-5
; Periodic boundary
pbc                 = xyz
DispCorr            = EnerPres
; Velocidades — regeneradas por réplica para independencia estadística
gen_vel             = yes
gen_temp            = ${TEMP}
gen_seed            = ${GEN_SEED}
EOF
    echo "  Production MDP written: ${PRODDIR}/prod.mdp (gen_seed=${GEN_SEED})"
else
    echo "  Using existing ${PRODDIR}/prod.mdp"
fi

# ── Preprocess (solo si el .tpr no existe aún) ──────────────────────────────
if [[ ! -f "${PRODDIR}/prod.tpr" ]]; then
    echo ""
    echo "  Preprocessing with $GMX grompp (rep ${REP})..."
    # Nota: NO se pasa -t equil/npt4.cpt — las velocidades se regeneran
    # (gen_vel=yes) a partir de las coordenadas equilibradas npt4.gro.
    $GMX grompp \
        -f "${PRODDIR}/prod.mdp" \
        -c equil/npt4.gro \
        -p "$GMX_TOP" \
        -n "$GMX_NDX" \
        -o "${PRODDIR}/prod.tpr" \
        -maxwarn 5 \
        2>&1 | tee "${PRODDIR}/prod_grompp.log"
    echo "  ✓ ${PRODDIR}/prod.tpr"
else
    echo "  ${PRODDIR}/prod.tpr already exists — skipping grompp."
fi

# ── Detect resume vs fresh start ────────────────────────────────────────────
CPI_FLAG=""
if [[ -f "${PRODDIR}/prod.cpt" ]]; then
    echo ""
    echo "  ► Detected ${PRODDIR}/prod.cpt — RESUMING from last checkpoint."
    CPI_FLAG="-cpi ${PRODDIR}/prod.cpt"
else
    echo ""
    echo "  ► No prod.cpt — FRESH START (réplica ${REP})."
    echo "  Target: ${PRODUCTION_NS} ns (${PRODUCTION_STEPS} steps)"
fi

# ── Run production ───────────────────────────────────────────────────────
# Banderas de GPU explícitas (antes se dejaba a la autodetección de GROMACS):
#   -nb gpu      : no-enlazante (Lennard-Jones/Coulomb corto alcance) en GPU
#   -pme gpu     : PME (electrostática larga distancia) en GPU
#   -bonded gpu  : términos enlazantes también en GPU
#   -gpu_id 0    : usa la GPU asignada por SLURM al job (ajusta si Yoltla
#                  expone más de una GPU por nodo/job)
# -update gpu (GPU-resident) puede dar otro salto de rendimiento en
# GROMACS 2021+ con LINCS + V-rescale + Parrinello-Rahman, pero pruébalo
# primero en una corrida corta — si el clúster/versión no lo soporta bien,
# quítalo y GROMACS seguirá corriendo con la actualización en CPU.
#
# -maxh 95.5: la partición "gpus" tiene un tope duro de --time=4-00:00:00
# (96 h). Con -maxh le decimos a GROMACS que pare LIMPIO y escriba
# checkpoint ~30 min antes de que SLURM mate el job a la fuerza, en vez de
# arriesgarnos a que lo corte a mitad de un paso. Si cambias --time arriba,
# ajusta -maxh a ese valor menos ~0.5 h.
echo ""
echo "  Running mdrun (réplica ${REP}) — GPU ..."
echo ""
srun -n 1 $GMX mdrun \
    -v \
    -deffnm "${PRODDIR}/prod" \
    -s "${PRODDIR}/prod.tpr" \
    -ntomp "$NCPU" \
    -nb gpu \
    -pme gpu \
    -bonded gpu \
    -gpu_id 0 \
    -cpt 30 \
    -maxh 95.5 \
    $CPI_FLAG \
    2>&1 | tee -a "${PRODDIR}/prod_mdrun.log"

# ── ¿Llegamos a los 100 ns, o hay que reencolar automáticamente? ───────────
# El tope de 4 días de la partición "gpus" casi seguro no alcanza para
# 100 ns en una sola corrida. En vez de que tengas que estar pendiente y
# resometer manualmente cada vez que se corte, el script se auto-reencola
# hasta que el checkpoint reporte haber llegado al paso objetivo.
CURRENT_STEP=$($GMX dump -cp "${PRODDIR}/prod.cpt" 2>/dev/null \
                | grep -m1 -E '^ *step ' | awk '{print $NF}')

if [[ -n "$CURRENT_STEP" && "$CURRENT_STEP" =~ ^[0-9]+$ && "$CURRENT_STEP" -ge "$PRODUCTION_STEPS" ]]; then
    echo ""
    echo "  ✓ Réplica ${REP}: simulación COMPLETA (step ${CURRENT_STEP}/${PRODUCTION_STEPS})."
else
    echo ""
    echo "  ► Réplica ${REP}: simulación INCOMPLETA (step actual: ${CURRENT_STEP:-desconocido} / ${PRODUCTION_STEPS})."
    echo "    Reencolando automáticamente para continuar desde el checkpoint..."
    sbatch \
        --job-name="prod100ns_rep${REP}" \
        --output="prod_rep${REP}_%j.out" \
        --error="prod_rep${REP}_%j.err" \
        "$0" "$REP"
fi

echo ""
echo "══════════════════════════════════════════════════"
echo "  Step 08 complete — Réplica ${REP}, Production MD finished!"
echo ""
echo "  Output files (en ${PRODDIR}/):"
echo "    prod.xtc     <- trajectory (compressed, 10 ps per frame)"
echo "    prod.edr     <- energy log"
echo "    prod.gro     <- final coordinates"
echo "    prod.cpt     <- checkpoint (also for resume)"
echo "    prod.log     <- MD log"
echo ""
echo "  Quick sanity checks:"
echo "    echo Temperature | $GMX energy -f ${PRODDIR}/prod.edr -o ${PRODDIR}/temp.xvg"
echo "    echo Density     | $GMX energy -f ${PRODDIR}/prod.edr -o ${PRODDIR}/dens.xvg"
echo "    $GMX rms -s ${PRODDIR}/prod.tpr -f ${PRODDIR}/prod.xtc -o ${PRODDIR}/rmsd.xvg"
echo ""
echo "  Suggested analyses (combinar las 3 réplicas):"
echo "    - RMSD de backbone de proteína (por réplica y promedio)"
echo "    - RMSD de glifosato y FAD"
echo "    - RMSF por residuo"
echo "    - Ocupación de puentes de H ($GMX hbond)"
echo "    - Energía libre de unión (g_mmpbsa — disponible en Yoltla)"
echo "══════════════════════════════════════════════════"