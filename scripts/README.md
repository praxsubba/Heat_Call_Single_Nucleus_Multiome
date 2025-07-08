# README.md `scripts`

## WGCNA Analysis Pipeline with DESeq2 VST Normalization and Cell-Type Deconvolution
This repository contains a comprehensive and well-annotated R pipeline for performing Weighted Gene Co-expression Network Analysis (WGCNA) on raw RNA-seq count data, normalized using DESeq2's Variance Stabilizing Transformation (VST). The pipeline includes gene length data processing, sample quality control, network construction, module detection, trait correlation, gene ontology enrichment, hub gene identification, network visualization, cell-type enrichment, and single-cell RNA-seq data deconvolution preparation. Additionally, a SLURM batch script for running CIBERSORTx deconvolution on a high-performance computing cluster is provided.

### Table of Contents
Setup and Data Preparation

Normalization and Quality Control

Network Construction and Module Detection

Module-Trait Relationships

Gene Ontology Enrichment

Hub Gene Identification and Validation

Network Visualization (Green Module)

Cell-Type Enrichment Analysis

Single-Cell RNA-seq Deconvolution Signature Preparation

Seurat Marker Gene Identification

SLURM Batch Script for CIBERSORTx Deconvolution


