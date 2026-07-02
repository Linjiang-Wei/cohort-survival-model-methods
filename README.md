# Data-free model methods

This repository contains data-free statistical model functions for time-to-event cohort analyses.

The scripts are intended for use only by researchers with permission to analyze the relevant data in an appropriate secure environment.

Only generic model functions are included.

The repository provides generic implementations of the following model-based steps:

1. Cox proportional hazards models for time-to-event outcomes.
2. Standardized risk-difference estimation by coefficient simulation.
3. Diagnosis-level outcome scans with false-discovery-rate correction.
4. Cross-fitted biomarker-profile adjustment.
5. Exposure-modifier interaction, joint-stratum and sensitivity analyses.

Users pass study-specific column names to the functions within their own secure environment.
