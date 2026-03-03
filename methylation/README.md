# Running Marmoset Workflow:

```
# Making Input jsons:
cd /private/groups/migalab/jmmenend/toil_runs/marmoset && \
    mkdir -p input_jsons && cd input_jsons

python3 /private/groups/migalab/jmmenend/toil_runs/scripts/launch_from_table.py \
    --data_table /private/groups/migalab/jmmenend/toil_runs/marmoset/mamoset_inputs.csv \
    --field_mapping /private/groups/migalab/jmmenend/toil_runs/marmoset/marmoset_input_mapping.csv \
    --workflow_name marmoset_workflow


# Launching Batch WDLs:
cd /private/groups/migalab/jmmenend/toil_runs/marmoset && \
    mkdir -p processing && cd processing && mkdir -p slurm_logs

sbatch \
    --job-name=marmoset_workflow \
    --array=[1-4]%2 \
    --partition=long \
    /private/groups/migalab/jmmenend/toil_runs/scripts/toil_sbatch_slurm.sh \
    --wdl /private/groups/migalab/jmmenend/toil_runs/wdls/marmoset_workflow.wdl \
    --sample_csv /private/groups/migalab/jmmenend/toil_runs/marmoset/mamoset_inputs.csv \
    --input_json_path '/private/groups/migalab/jmmenend/toil_runs/marmoset/input_jsons/${SAMPLE_ID}_marmoset_workflow.json'


# Make Output JSON
cd /private/groups/migalab/jmmenend/toil_runs/marmoset

python3 /private/groups/migalab/jmmenend/toil_runs/scripts/update_table_with_outputs.py \
    --input_data_table /private/groups/migalab/jmmenend/toil_runs/marmoset/mamoset_inputs.csv \
    --output_data_table /private/groups/migalab/jmmenend/toil_runs/marmoset/mamoset_outputs.csv \
    --json_location '/private/groups/migalab/jmmenend/toil_runs/marmoset/processing/{sample_id}/{sample_id}_marmoset_workflow_outputs.json'
```


