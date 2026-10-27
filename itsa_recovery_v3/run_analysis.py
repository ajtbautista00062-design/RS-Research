"""Explicit V3 entry point. Importing this module never runs the pipeline."""
from pathlib import Path
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(1,str(HERE/"preserved_v1"))
import argparse,hashlib,importlib.util,json,logging,shutil,unittest
from datetime import datetime,timezone
import config as cfg
from sarflow.temporal import discover,local_path,file_digest
from model import design
from validation import validate_settings

def synthetic_validation():
    suite=unittest.defaultTestLoader.discover(str(HERE/"tests"),pattern="test_*.py")
    result=unittest.TextTestRunner(verbosity=1).run(suite)
    if not result.wasSuccessful():
        raise RuntimeError("Synthetic/software validation failed; no data processing allowed.")
    return result.testsRun

def load_legacy_config():
    path=HERE/"preserved_v1"/"config.py"
    spec=importlib.util.spec_from_file_location("v1_settings_snapshot",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

def stage_gate(stage,obs,pilot):
    counts=tuple(sum(o.phase==p for o in obs) for p in ("pre","event","post"))
    expected=cfg.EXPECTED_PILOT_COUNTS if stage=="pilot" else cfg.EXPECTED_FULL_COUNTS
    if counts!=expected:
        raise ValueError(f"Stage {stage} requires {expected} unique pairs; found {counts}.")
    if cfg.ZERO_SEMANTICS not in ("unresolved","valid","background"):
        raise ValueError("Invalid ZERO_SEMANTICS.")
    if cfg.ZERO_SEMANTICS!="unresolved" and not cfg.ZERO_SEMANTICS_JUSTIFICATION.strip():
        raise ValueError("Document the source evidence before declaring zero semantics resolved.")
    if stage=="full":
        if not cfg.ALLOW_FULL_DATASET or not pilot:
            raise ValueError("Full execution requires ALLOW_FULL_DATASET and --validated-pilot.")
        folder=local_path(pilot,cfg.PROJECT_ROOT)
        state=json.loads((folder/"run_state.json").read_text())
        manifest=json.loads((folder/"run_manifest.json").read_text())
        if (state.get("status")!="complete" or not manifest.get("technical_validation_passed")
                or manifest.get("stage")!="pilot" or manifest.get("implementation_version")!="3"):
            raise ValueError("A technically validated completed V3 pilot is required.")
        if cfg.ZERO_SEMANTICS=="unresolved" or not cfg.ACQUISITION_COMPARABILITY_CONFIRMED:
            raise ValueError("Resolve zero semantics and acquisition comparability before full execution.")
    return counts

def provenance(obs,stage,tests,counts):
    settings=vars(cfg.analysis_settings())
    sources={}
    for path in HERE.rglob("*.py"):
        if "__pycache__" not in path.parts:
            sources[str(path.relative_to(HERE))]=hashlib.sha256(path.read_bytes()).hexdigest()
    rows=[dict(name=o.name,path=str(o.path),phase=o.phase,start=o.start.isoformat(),end=o.end.isoformat(),
               reference=o.reference.isoformat(),included_primary_itsa=o.phase!="event",sha256=file_digest(o.path)) for o in obs]
    return dict(implementation_version="3",stage=stage,counts=counts,synthetic_tests_passed=tests,observations=rows,
        boundary_sha256=file_digest(Path(cfg.Study_Area)),
        zero_semantics=cfg.ZERO_SEMANTICS,zero_semantics_justification=cfg.ZERO_SEMANTICS_JUSTIFICATION,
        comparability_confirmed=cfg.ACQUISITION_COMPARABILITY_CONFIRMED,
        analysis_status="EXPLORATORY - operational thresholds and error assumptions require validation",
        authoritative=False,primary_spatial_fill=False,event_origin=cfg.EVENT_ORIGIN,
        reference_date=cfg.REFERENCE_DATE,thresholds=settings,source_code_sha256=sources,
        classification="ITSA trajectory-derived with observed confirmation and event evidence",
        uncertainty="Operational guards; not calibrated joint confidence intervals",
        comparisons="Existing V1 raster may reflect different masking/correction; transitions are descriptive")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--stage",choices=("pilot","full"),default="pilot")
    parser.add_argument("--validated-pilot",type=Path)
    args=parser.parse_args()
    tests=synthetic_validation()
    logging.basicConfig(level=logging.INFO,format="%(asctime)s %(levelname)s %(message)s")
    logging.info("Discovering distinct real inputs; no synthetic coherence is used in analysis.")
    obs,duplicates=discover(cfg)
    counts=stage_gate(args.stage,obs,args.validated_pilot)
    calibration=None
    if cfg.RECALIBRATE_FROM_PRE:
        from calibration import calibrate
        logging.info("PRE-only unfilled noise-calibration pass; no class-frequency optimization.")
        calibration=calibrate(cfg,obs)
        cfg.CALIBRATED_LOSS_MIN=calibration["operational_floor"]
    validate_settings(cfg.analysis_settings())
    design(obs,cfg.EVENT_ORIGIN)
    dates=f"{min(o.start for o in obs):%Y%m%d}_{max(o.end for o in obs):%Y%m%d}"
    stamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root=local_path(cfg.Output_Folder,cfg.PROJECT_ROOT)/f"itsa_v3_{args.stage}_{dates}_{stamp}"
    root.mkdir(parents=True,exist_ok=False)
    log=logging.FileHandler(root/"processing.log",encoding="utf-8")
    log.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logging.getLogger().addHandler(log)
    state=root/"run_state.json"
    state.write_text(json.dumps(dict(status="running")),encoding="utf-8")
    metadata=provenance(obs,args.stage,tests,counts)
    metadata["loss_calibration"]=calibration
    try:
        if calibration:
            (root/"loss_calibration.json").write_text(json.dumps(calibration,indent=2),encoding="utf-8")
        for name in ("config.py","METHOD.md","TECHNICAL_REPORT.md","original_source_hashes.json"):
            shutil.copy2(HERE/name,root/name)
        from workflow import execute
        tags=dict(analysis_status=metadata["analysis_status"],zero_semantics=cfg.ZERO_SEMANTICS,
                  source_boundary=cfg.Study_Area,event_origin=cfg.EVENT_ORIGIN,pair_reference=cfg.REFERENCE_DATE,
                  thresholds=json.dumps(metadata["thresholds"]),primary_correction="none",
                  stage=args.stage,run_id=root.name)
        details=execute(cfg,obs,root,tags,load_legacy_config())
        metadata.update(details,duplicate_files=duplicates,
                        technical_validation_passed=details["independent_real_checks"]>0)
        (root/"run_manifest.json").write_text(json.dumps(metadata,indent=2,default=str),encoding="utf-8")
        state.write_text(json.dumps(dict(status="complete",technical_validation_passed=metadata["technical_validation_passed"])),
                         encoding="utf-8")
        logging.info("Completed exploratory outputs: %s",root)
    except (Exception,KeyboardInterrupt) as exc:
        logging.exception("Run stopped; partial artifacts are not validated outputs.")
        state.write_text(json.dumps(dict(status="interrupted" if isinstance(exc,KeyboardInterrupt) else "failed",
                                         error_type=type(exc).__name__,error=str(exc))),encoding="utf-8")
        raise

if __name__=="__main__":
    main()
