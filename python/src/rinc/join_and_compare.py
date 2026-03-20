'''
Created on Jan 5, 2026
Gather gap, hgvs, and annovar information; perform comparison, and write results to csv.
  
@author: pleyte
'''
import argparse
import logging.config
from rinc.util.log_config import LogConfig
import pandas as pd
from functools import reduce
from rinc.variant_transcript import VariantTranscript
from enum import Enum
import itertools
import numpy as np

class NomenclatureTools(Enum):
    TFX = "tfx"
    ANNOVAR = "annovar"
    SNPEFF = "snpeff"
    CGD = "cgd"
    VEP_REFSEQ = "vep_refseq"
    VEP_HG19 = "vep_hg19"
    VARIANT_VALIDATOR = "vv"
    MUTALYZER = "mutalyzer"
    REFERENCE_GAP = "gap" # Not actually a nomenclature tool
    PREFERRED_TRANSCRIPT = "cgd_preferred_transcript" # Not actually a nomenclature tool

class JoinAndCompare(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._index_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
                
    def _get_all_variant_transcripts(self, dataframes: list[pd.DataFrame]) -> list[VariantTranscript]:
        """
        Return a list of every variant transcript in all of the dataframes 
        """
        join_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        
        merged_df = reduce(lambda left, right: pd.merge(left, right, on=join_cols, how='outer'), dataframes)
        
        variant_transcripts = []
        for row in merged_df[['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']].itertuples(index=False):
            variant_transcripts.append(VariantTranscript(row.chromosome, row.position, row.reference, row.alt, row.cdna_transcript))

        return variant_transcripts
    
    def get_info_df(self, csv_file: str, info_tool_name: str):
        """
        Open and return a csv dataframe adding making it a MultiIndex with 0 level id of info_tool_name        
        """
        if not csv_file:
            return None
        
        df = pd.read_csv(csv_file)
        df.attrs['info_tool'] = info_tool_name
        return df
        
    def get_nomenclature_df(self, nomenclature_tool_name: str, csv_file: str, remove_label: str=None):
        """
        Read a file containing variant transcript nomenclature.
        If the file doesn't exist then an empty dataframe is returned.
        
        The dataframe columns all have identifiers in the (eg tfx.exon or p_dot1.annovar) but they get in the way 
        here so the `removelabel` string is removed from the columns. 
        """
        df = pd.read_csv(csv_file, dtype=str)

        index_cols = self._index_cols
        df.set_index(index_cols, inplace=True, drop=False)

        if remove_label:
            self._logger.info(f"Removing {remove_label} label from fields in {nomenclature_tool_name} df")
            df.columns = [x.replace(remove_label, '') for x in df.columns]
            
        df.columns = pd.MultiIndex.from_product([[nomenclature_tool_name], df.columns])
        
        df.attrs['nomenclature_tool'] = nomenclature_tool_name
        self._logger.info(f"Read {df.shape[0]} rows for {nomenclature_tool_name} from {csv_file}")
        return df

    def _get_variant_transcript_row(self, df, chromosome, position, reference, alt, cdna_transcript):
        """
        Query a dataframe for a row with the variant and transcript
        """
        create_mask = lambda df, c, p, r, a, t: (
            (df['chromosome'] == c) & 
            (df['position'] == p) & 
            (df['reference'] == r) & 
            (df['alt'] == a) &
            (df['cdna_transcript'] == t))
        
        df_row = df[create_mask(df, chromosome, position, reference, alt, cdna_transcript)]
        return df_row
    
    def _get_merged_df(self, dataframes: list[pd.DataFrame]):
        """
        Merge all the dataframes into one. 
        I don't want to inner join them all because that is too limiting.
        Outter join ends up with a lot of sparse fields, but it's what i've decided to do.
        """        
        merged_df = reduce(
            lambda left, right: pd.merge(left, right, on=['chromosome','position','reference', 'alt','cdna_transcript'], how='outer'), 
            dataframes
        )

        return merged_df    
        
    def _remove_rows_without_cdna_transcript_version(self, merged_df: pd.DataFrame):
        """
        Remove rows where the cdna transcript doesn't have a version. CGD should be the only 
        datasource that has transcripts w/o accession version. 
        """        
        df_filtered = merged_df[merged_df.index.get_level_values('cdna_transcript').str.contains('.', regex=False)]
        self._logger.info(f"Removed {len(merged_df) - len(df_filtered)} rows where transcript accession does not have version")
        return df_filtered
    
    def _remove_rows_without_cgd_and_tfx_datasource_values(self, merged_df: pd.DataFrame, dataframes: list):
        """
        Our goal is to validate transcript effects and cgd data. If neither cgd or tfx have informatio for a particular 
        transcript then we have no use for it. 
        """
        cgd_df = next(x_df for x_df in dataframes if x_df.attrs.get('nomenclature_tool') == NomenclatureTools.CGD.value)
        tfx_df = next(x_df for x_df in dataframes if x_df.attrs.get('nomenclature_tool') == NomenclatureTools.TFX.value)
        assert merged_df.index.names == cgd_df.index.names == tfx_df.index.names, "all three dataframes are expected to have the same index"
                
        cgd_index = cgd_df.index
        tfx_index = tfx_df.index

        # Combine them (Union). This creates a single index containing every variant from both sources
        # Using union() is safer than concat because it handles MultiIndex levels correctly
        allowed_index = cgd_index.union(tfx_index)
        
        # Filter the merged_df
        final_df = merged_df[merged_df.index.isin(allowed_index)]
        
        self._logger.info(f"Removed {len(merged_df) - len(final_df)} rows that have neither cgd or tfx values")
        return final_df
    
    def _fill_in_missing_genomic_variant_id(self, df):
        '''
        The variant transcripts from CGD have a genomic_variant_id. Transcripts that aren't known to cgd will have a blank genomic_variant_id.
        This function fills in the blank genomic_variant_id.
        '''
        if ('cgd', 'genomic_variant_id') not in df.columns:
            raise ValueError("cgd and genomic_variant_id are not in this dataframe")

        gv_id_col = ('cgd', 'genomic_variant_id')
        variant_levels = ['chromosome', 'position', 'reference', 'alt']
        
        # Group by the variant levels and use transform to propagate the ID
        # ffill (forward fill) and bfill (backward fill) inside the group 
        # ensures the ID spreads to all rows regardless of where it started.
        df[gv_id_col] = (df.groupby(level=variant_levels)[[gv_id_col]].transform(lambda x: x.ffill().bfill()))
        
        # Optional: Verify how many IDs were successfully recovered
        total_rows = len(df)
        remaining_nan = df[gv_id_col].isna().sum()
        self._logger.info(f"Genomic variant id fill complete. {total_rows - remaining_nan} rows now have an ID. {remaining_nan} rows still missing an ID (no match found in CGD)")
        return df
        
    def get_comparison_df(self, dataframes: list[pd.DataFrame], gff_and_uta_exon_gap_info_df, preferred_transcript_df=None):
        """
        Compare dataframes with each other and return a dataframe with comparison score         
        """
        merged_df = self._get_merged_df(dataframes)
                            
        # Remove unwanted rows 
        merged_df = self._remove_rows_without_cgd_and_tfx_datasource_values(merged_df, dataframes)        
        merged_df = self._remove_rows_without_cdna_transcript_version(merged_df)
        
        merged_df = self._calculate_pairwise_score(merged_df, dataframes, fields_to_compare=['g_dot'], new_field_name='g_dot_concordance')
        merged_df = self._calculate_pairwise_score(merged_df, dataframes, fields_to_compare=['c_dot'], new_field_name='c_dot_concordance')
        merged_df = self._calculate_pairwise_score(merged_df, dataframes, fields_to_compare=['p_dot1'], new_field_name='p_dot_cordance')
        merged_df = self._calculate_pairwise_score(merged_df, dataframes, fields_to_compare=['c_dot', 'p_dot1'], new_field_name='c+p_concordance')
        merged_df = self._calculate_pairwise_score(merged_df, dataframes, fields_to_compare=['exon', 'c_dot', 'p_dot1'], new_field_name='c+p+exon_concordance')
        merged_df = self._add_vep_refseq_ref_mismatch_field(merged_df, new_field_name='vep_refseq_mismatch')
        
        # Fill in teh genomic_variant_id field 
        merged_df = self._fill_in_missing_genomic_variant_id(merged_df)
        
        # Add a field indicating when cgd & annovar agree but disagree with tfx & vepHg19 on c_dot
        merged_df = self._add_consensus_conflict_field(merged_df, ['cgd', 'annovar'], ['tfx', 'vep_hg19'], ['c_dot'], 'ca_vs_tv_conflict')
        
        # Add a field indicating which transcripts are known to have gaps/misalignment
        merged_df = self._add_gap_and_cigar_field(merged_df, gff_and_uta_exon_gap_info_df)

        # Add a cgdAnnovar_vs_tfx field
        merged_df = self._add_cgd_and_annovar_vs_tfx_field(merged_df)
        
        # Add a is_variant_tx_gap
        merged_df = self._add_variant_tx_gap_field(merged_df)
        
        # Add a field inndicating which accessions are preferred/reported transcripts
        if preferred_transcript_df is not None:
            merged_df = self._add_preferred_transcript_field(merged_df, preferred_transcript_df)
                
        return merged_df

    def _add_variant_tx_gap_field(self, df: pd.DataFrame):
        """
        Add a field that indicates a variant is associated with a transcript that has a reference-refseq mismatch.
        This field will be true for all variant-transcripts where the variant is associated with a reference-refseq mismatch even if the cdna_transcript is not the transcript that has the mismatch.  
        """
        df[('gap', 'is_variant_tx_gap')] = (
            df[('gap', 'is_gap_transcript')]
            .groupby(level=['chromosome', 'position', 'reference', 'alt'])
            .transform('max')
        )
        
        original_gaps = df[('gap', 'is_gap_transcript')].sum()
        variant_gaps = df[('gap', 'is_variant_tx_gap')].sum()

        self._logger.info(f"Transcript gaps: {original_gaps}")
        self._logger.info(f"Total rows marked by variant gap: {variant_gaps}")
        
        return df
        
    def _get_concordance_counts(self, df: pd.DataFrame, fields: list):
        """
        Generate concordage table
        """
        sources = [NomenclatureTools.CGD.value, NomenclatureTools.TFX.value] # NomenclatureTools.ANNOVAR.value,

        # 1. Build a mask where all sources have valid data for ALL metrics provided
        # This ensures we aren't comparing a 'c_dot' if the 'p_dot' is missing
        presence_masks = []
        for m in fields:
            m_mask = (df[sources].xs(m, axis=1, level=1).notna() & 
                      (df[sources].xs(m, axis=1, level=1) != '')).all(axis=1)
            presence_masks.append(m_mask)
        
        # Combine all presence masks (must have data for every metric)
        final_mask = pd.concat(presence_masks, axis=1).all(axis=1)
        df_filtered = df[final_mask]
        
        # 2. Build the count matrix
        matrix = pd.DataFrame(index=sources, columns=sources)
        
        for s1 in sources:
            for s2 in sources:
                # Check for equality across ALL metrics in the list
                # .all(axis=1) ensures row matches only if c_dot AND p_dot (etc) match
                matches = (df_filtered[s1][fields] == df_filtered[s2][fields]).all(axis=1).sum()
                matrix.loc[s1, s2] = matches
        
        index_df = df_filtered.index.to_frame(index=False)
        variant_coords = ['chromosome', 'position', 'reference', 'alt']
        distinct_variant_count = len(index_df[variant_coords].drop_duplicates())
        
        # Return variant count, transcript count, and the matrix 
        return distinct_variant_count, len(df_filtered), matrix.astype(int)
        
    def _generate_concordance_tables(self, df: pd.DataFrame):
        """
        Generate concordage tables (aka Co-occurrence Matrices) for CGD and Tfx agreement on c., p., and both
        """
        
        variant_count_c, transcripts_count_c, table_c = self._get_concordance_counts(df, ['c_dot'])
        print("-" * 50)
        print(f"c. concordance: Analysis based on {variant_count_c} variants and {transcripts_count_c} transcripts where all sources have c dot.")        
        print(table_c)
        print(f"")
        
        variant_count_p, transcripts_count_p, table_p = self._get_concordance_counts(df, ['p_dot1'])
        print("-" * 50)
        print(f"p. concordance: Analysis based on {variant_count_p} variants and {transcripts_count_p} transcripts where all sources have p dot.")
        print(table_p)
        print(f"")
        
        variant_count_cp, transcripts_count_cp, table_cp = self._get_concordance_counts(df, ['c_dot', 'p_dot1'])
        print("-" * 50)
        print(f"c. and p. concordance: Analysis based on {variant_count_cp} variants and {transcripts_count_cp} transcripts where all sources have c dot and p dot.")
        print(table_cp)
        print(f"")
        
    def _add_preferred_transcript_field(self, merged_df, preferred_transcript_df):
        """
        Add a field with a "1" wherever the cdna_transcript matches in mergrd_df matches a cdna_transcrpt in preferred_tranascript_df 
        """
        pt_label = NomenclatureTools.CGD.value

        # Use .isin() - this is a vectorized membership check
        # It checks all rows at once in optimized C code
        is_pref_bool = merged_df.index.get_level_values('cdna_transcript').isin(
            preferred_transcript_df['cdna_transcript']
        )
        
        merged_df[(pt_label, 'is_preferred')] = is_pref_bool.astype(int)
        return merged_df
        
    def _add_gap_and_cigar_field(self, merged_df, gff_and_uta_exon_gap_info_df):
        """
        Join the dataframe that is a list of all the transcripts which are known to have reference gaps.
        The dataframe has two cigar columns: gff_cigars, uta_cigars
        Add a field is_gap_transcript to make it clear which transcripts were joined.
        """
        gap_label = NomenclatureTools.REFERENCE_GAP.value
        
        # 1. Prepare the lookup table: set 'accession' as the index
        # Ensure it's unique so the mapping is 1:1
        lookup = gff_and_uta_exon_gap_info_df.drop_duplicates('accession').set_index('accession')
    
        # 2. Extract the 'cdna_transcript' level from the merged_df index
        # This gets the transcript ID for every row without moving it out of the index
        transcript_ids = merged_df.index.get_level_values('cdna_transcript')
    
        # 3. Map the cigar values directly to new MultiIndex columns
        # .map() follows the index level values and pulls data from the lookup table
        merged_df[(gap_label, 'gff_cigars')] = transcript_ids.map(lookup['gff_cigars'])
        merged_df[(gap_label, 'uta_cigars')] = transcript_ids.map(lookup['uta_cigars'])
    
        # 4. Create your binary flag
        has_gff = merged_df[(gap_label, 'gff_cigars')].notna()
        has_uta = merged_df[(gap_label, 'uta_cigars')].notna()
        merged_df[(gap_label, 'is_gap_transcript')] = (has_gff | has_uta).astype(int)
    
        return merged_df
        
    def _add_cgd_and_annovar_vs_tfx_field(self, df: pd.DataFrame):
        """
        Add a field that indicates which rows have the same on c. and p. for CGD and Annovar but Tfx has a different c. or p.  
        """
        sources = [NomenclatureTools.CGD.value, NomenclatureTools.ANNOVAR.value, NomenclatureTools.TFX.value]
        
        # 1. Create a mask where all three sources have non-empty/non-null values
        presence_masks = [
            (df[(s, 'c_dot')].notna() & (df[(s, 'c_dot')] != '')) & 
            (df[(s, 'p_dot1')].notna() & (df[(s, 'p_dot1')] != ''))
            for s in sources
        ]
        all_present = reduce(lambda x, y: x & y, presence_masks)
        
        # 2. Logic: CGD cdot and pdot  matches Annovar
        cgd_and_annovar_match = (
            (df[(NomenclatureTools.CGD.value, 'c_dot')] == df[(NomenclatureTools.ANNOVAR.value, 'c_dot')]) & 
            (df[(NomenclatureTools.CGD.value, 'p_dot1')] == df[(NomenclatureTools.ANNOVAR.value, 'p_dot1')])
        )
        
        # 3. Logic: TFX matches CGD (and by extension Annovar if cgd_ann_match is True)
        tfx_match = (
            (df[(NomenclatureTools.TFX.value, 'c_dot')] == df[(NomenclatureTools.CGD.value, 'c_dot')]) & 
            (df[(NomenclatureTools.TFX.value, 'p_dot1')] == df[(NomenclatureTools.CGD.value, 'p_dot1')])
        )
        
        # 4. Define Conditions
        # All match = all present AND (cgd==annovar) AND (tfx matches)
        cond_all_match = all_present & cgd_and_annovar_match & tfx_match

        # TFX differs = all present AND (cgd==annovar) AND (tfx doesn't match)
        cond_tfx_outlier = all_present & cgd_and_annovar_match & ~tfx_match
        
        # 5. Assign to the multi-index
        df[('scores', 'cgdAnnovar_vs_tfx')] = np.select(
            [cond_all_match, cond_tfx_outlier], 
            [1, 0], 
            default=-1
        )
        
        return df
        
        
    def _calculate_pairwise_score(self, merged_df: pd.DataFrame, tool_dataframes: list[pd.DataFrame], fields_to_compare: list[str], new_field_name):
        """
        Compare the fields_to_compare   
        """
        tool_dataframe_names = [x.attrs['nomenclature_tool'] for x in tool_dataframes]
        
        pairs = list(itertools.combinations(tool_dataframe_names, 2))
        
        total_matches = pd.Series(0, index=merged_df.index)
                
        possible_comparisons = pd.Series(0, index=merged_df.index)
        
        for t1, t2 in pairs:
            
            t1_has_fields = all((t1, f) in merged_df.columns for f in fields_to_compare)
            t2_has_fields = all((t2, f) in merged_df.columns for f in fields_to_compare)
            
            # If the datasource doesn't provide a field we can't compare it (eg CGD does not have a g_dot)
            if not t1_has_fields:                
                self._logger.info(f"Unable to calculate pairwise score of {fields_to_compare} because {t1} does not have that one of the fields")                
                continue
            elif not t2_has_fields:
                self._logger.info(f"Unable to calculate pairwise score of {fields_to_compare} because {t2} does not have that one of the fields")                
                continue
                 
            # 1. Determine if BOTH tools have data for ALL fields in the list
            # We start with True and 'AND' it with each field's presence
            both_have_all_data = pd.Series(True, index=merged_df.index)
            for field in fields_to_compare:
                both_have_all_data &= (merged_df[t1, field].notna() & merged_df[t2, field].notna())
            
            # 2. Determine if ALL fields match exactly between the two tools
            # We start with True and 'AND' it with each field's comparison
            all_fields_match = pd.Series(True, index=merged_df.index)
            for field in fields_to_compare:
                all_fields_match &= (merged_df[t1, field] == merged_df[t2, field])
            
            # 3. A 'Point' is awarded only if they match AND both had data
            match_mask = all_fields_match & both_have_all_data
            
            # 3. Add to counters
            total_matches += match_mask.astype(int)
            possible_comparisons += both_have_all_data.astype(int)
            
        # 4. Calculate Score: matches / comparisons 
        # (avoid division by zero with .replace)        
        merged_df['scores', new_field_name] = total_matches / possible_comparisons.replace(0, np.nan)
        return merged_df
    
        
    def _get_field_name(self, pattern, df: pd.DataFrame, optional=False):
        """
        Return the name of the field that has c. in it
        """
        field_names = []
        for x in df.columns:
            if pattern in x:
                field_names.append(x)
    
        if len(field_names) > 1: 
            raise ValueError(f"Multiple field names matching {pattern}: {field_names}")
        if not optional and not field_names:
            raise ValueError(f"Unable to find field matching {pattern} in {df.attrs['nomenclature_tool']}")
        elif field_names:
            return field_names[0]
        else:
            return None
    
    def _add_vep_refseq_ref_mismatch_field(self, merged_df: pd.DataFrame, new_field_name: str):
        """
        Add a field that is 1 when vepRefSeq GIVEN_REF != USED_REF
        """
        mismatch_mask = (
            (merged_df[(NomenclatureTools.VEP_REFSEQ.value, 'GIVEN_REF')] != merged_df[(NomenclatureTools.VEP_REFSEQ.value, 'USED_REF')]) & 
            merged_df[(NomenclatureTools.VEP_REFSEQ.value, 'GIVEN_REF')].notna() & 
            merged_df[(NomenclatureTools.VEP_REFSEQ.value, 'USED_REF')].notna()
            )
        
        merged_df[('scores', new_field_name)] = mismatch_mask.astype(int)
        return merged_df
        
    def _get_row(self, nomenclature_tool, rows: list):
        """
        Return a row where the  nomenclature_tool attribute matchces the parameter  
        """
        for x in rows:
            if x.attrs['nomenclature_tool'] == nomenclature_tool:
                return x
        
        return None
        
    def _write_sheet_raw(self, writer, workbook, df, sheet_name):
        """
        Add the raw data worksheet 
        """
        df.to_excel(writer, sheet_name=sheet_name)
        worksheet = writer.sheets[sheet_name]
        
        bold_fmt = workbook.add_format({'bold': True, 'bg_color': '#D3D3D3', 'border': 1})
        worksheet.set_row(0, None, bold_fmt)
        worksheet.set_row(1, None, bold_fmt)
        
        worksheet.freeze_panes(2, 0)
        worksheet.set_column('A:ZZ', 18)

    def _write_sheet_comparison(self, writer, workbook, flat_df: pd.DataFrame, sheet_name):
        """
        Format the data to make it easier to read  
        """
        flat_df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # Bold the single header row in the new sheet
        summary_sheet = writer.sheets[sheet_name]
        
        bold_fmt = workbook.add_format({'bold': True, 'bg_color': '#D3D3D3', 'border': 1})
        summary_sheet.set_row(0, None, bold_fmt)
        
        summary_sheet.freeze_panes(1, 0)
        summary_sheet.set_column('A:E', 12) # Genomic coords
        summary_sheet.set_column('F:ZZ', 25) # Nomenclature fields
        
    def _write_sheet_summary(self, writer, workbook, flat_df, sheet_name):
        """
        """
        # 1. Identify all your score columns
        score_cols = [c for c in flat_df.columns if 'scores_' in c]
        
        summary_data = []
        
        for col in score_cols:
            # Get the data for this specific score
            data = flat_df[col]
            
            # How many variants had enough data to actually produce a score?
            total_comparable = data.notna().sum()
            
            # How many of those were perfect matches (score == 1.0)?
            perfect_matches = (data == 1.0).sum()
            
            # Calculate percentage
            percentage = (perfect_matches / total_comparable * 100) if total_comparable > 0 else 0
            
            summary_data.append({
                'Metric': col.replace('scores_', ''),
                'Perfect Concordance (n)': perfect_matches,
                'Total Comparable Variants': total_comparable,
                'Concordance Rate (%)': round(percentage, 2)
            })
        
        # 2. Create the Summary DataFrame
        metrics_df = pd.DataFrame(summary_data)
        
        # 3. Write to a new sheet in your Excel file
        # (Add this inside your 'with pd.ExcelWriter...' block)
        metrics_df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # 4. Optional: Add a little formatting to the Metrics sheet
        ws_metrics = writer.sheets[sheet_name]
        ws_metrics.set_column('A:A', 35) # Make the Metric name wide
        ws_metrics.set_column('B:D', 20) # Center the numbers
        
    def _get_flattened_data_frame(self, df):
        """
        """
        # We reset the index to turn the 5 genomic levels (Chr, Pos, etc.) into regular columns        
        flat_df = df.copy().reset_index()
        
        # 2. Force everything to string, handling MultiIndex tuples specifically
        new_cols = []
        for c in flat_df.columns:
            if isinstance(c, tuple):
                # Join the tuple parts with underscore, filter out empty strings
                clean_col = "_".join([str(part) for part in c if str(part).strip()])
                new_cols.append(clean_col)
            else:
                new_cols.append(str(c))
        
        
        flat_df.columns = new_cols
        
        # Identify columns index columns
        index_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
    
        # every dataframe brought in its own index column (eg cgd_chromosome) and we don't need to see those repeated.  
        redundant_suffixes = [
            'chromosome', 'position', 'reference', 'alt', 'cdna_transcript'
        ] 
           
        # nom_cols are all the not index "nomenclature" and other fields 
        nom_cols = [c for c in flat_df.columns 
                    if c not in index_cols
                    and not any(c.endswith(suffix) for suffix in redundant_suffixes)
                    ]
                                
        field_order = ['c_dot', 'p_dot1', 'exon', 'g_dot']
        summary_flags = ['cgd_is_preferred', 'gap_is_gap_transcript']

        # 5. Perform the Tiered Sort
        # Tier 0: nomenclature ending in field_order
        # Tier 1: summary flags
        # Tier 2: concordance scores
        # Tier 3: everything else
        sorted_nom_cols = sorted(
            nom_cols, 
            key=lambda x: (
                0 if any(x.endswith(f) for f in field_order) else
                1 if x in summary_flags else
                2 if x.startswith('scores_') else                
                3,
                # Sub-sort Tier 0 based on field_order index
                next((i for i, f in enumerate(field_order) if x.endswith(f)), 0),
                # Final alphabetical tie-breaker,
                x
            )
        )        
        
        final_column_order = index_cols + sorted_nom_cols
        return flat_df[final_column_order]
        
    def _add_consensus_conflict_field(self, merged_df: pd.DataFrame, list_a: list[str], list_b: list[str], fields: list[str], new_field_name):
        """
        Identifies rows where ListA tools agree, ListB tools agree, but Group A disagrees with Group B.
        """
        # Initialize a mask that starts as True for all rows
        # We will 'AND' this with our conditions
        final_mask = pd.Series(True, index=merged_df.index)
    
        for field in fields:
            # --- 1. Check Internal Consensus for List A ---
            # All must be non-null and equal to the first tool in the list
            a_first = list_a[0]
            a_consensus = pd.Series(True, index=merged_df.index)
            
            for tool in list_a:
                # Must have data
                a_consensus &= merged_df[tool, field].notna()
                # Must match the first tool
                a_consensus &= (merged_df[tool, field] == merged_df[a_first, field])
            
            # --- 2. Check Internal Consensus for List B ---
            b_first = list_b[0]
            b_consensus = pd.Series(True, index=merged_df.index)
            
            for tool in list_b:
                b_consensus &= merged_df[tool, field].notna()
                b_consensus &= (merged_df[tool, field] == merged_df[b_first, field])
                
            # --- 3. Check for Disagreement between Group A and Group B ---
            groups_disagree = (merged_df[a_first, field] != merged_df[b_first, field])
            
            # Combine for this specific field
            field_criteria = a_consensus & b_consensus & groups_disagree
            
            # Update the final mask (Row must meet criteria for ALL fields, e.g., c_dot AND p_dot)
            final_mask &= field_criteria
    
        self._logger.info(f"Consensus conflict for {list_a} vs {list_b} has {final_mask.sum()} rows")
        
        # Add a boolean field indicating rows that match the criteria
        merged_df[('scores', new_field_name)] = final_mask
        return merged_df
        

    def write(self, 
              out_file, 
              comparison_df: pd.DataFrame, 
              include_raw_sheet: bool = False):
        """
        Write the dataframe with its new comparison fields to file
        #merged_df.columns = [f"{source}_{field}" for source, field in merged_df.columns]
        #merged_df.to_csv("/tmp/merged_pairwise_df.csv", index=True)  
        """                        
        #sort_cols = ['c+p+exon_concordance', 'chromosome', 'position', 'reference', 'alt']
        #sort_orders = [False, True, True, True, True]
        comparison_df.sort_values(
            by=('scores', 'c+p+exon_concordance'), 
            ascending=False, 
            inplace=True
            )
        
        self._logger.info("Creating workbook")
        
        with pd.ExcelWriter(out_file, engine='xlsxwriter') as writer:
            workbook  = writer.book
            
            flatened_df = self._get_flattened_data_frame(comparison_df)
            
            self._write_sheet_summary(writer, workbook, flatened_df, 'Summary')            
            self._write_sheet_comparison(writer, workbook, flatened_df, 'Comparison')
            
            if include_raw_sheet:
                self._write_sheet_raw(writer, workbook, comparison_df, 'Raw')
            

        self._logger.info(f"Wrote data frame with {comparison_df.shape[0]} rows and {comparison_df.shape[1]} columns to {out_file}")
        
def _parse_args():
    parser = argparse.ArgumentParser(description='Read Annovar multianno file, extract values and write to new csv')

    parser.add_argument('--tfx_nomenclature', help="Optional Transcript Effect/tfx nomenclature (csv)", required=False)
    parser.add_argument('--cgd_nomenclature', help='CGD nomenclature (csv)', required=False)
    parser.add_argument('--variant_validator_nomenclature', help='Variant Validator nomenclature (csv)', required=False)
    parser.add_argument('--mutalyzer_nomenclature', help='Mutalyzer nomenclature (csv)', required=False)
    parser.add_argument('--snpeff_nomenclature', help='SnpEff nomenclature (csv)', required=False)
    
    parser.add_argument('--annovar_nomenclature', help='Annovar nomenclature (csv)', required=True)    
    parser.add_argument('--vep_refseq_nomenclautre', help='VEP Refseq nomenclature (csv)', required=True)
    parser.add_argument('--vep_hg19_nomenclature', help='Vep hg19 nomenclature (csv)', required=True)

    parser.add_argument('--gff_and_uta_exon_gap_info', help='Transcripts with refseq/hg19 alignment gaps (csv)', required=True)
    parser.add_argument('--preferred_transcripts', help='Preferred/reported transcripts (csv)', required=False)
    
    parser.add_argument("--out", help="output file (xlsx)", required=True)
    
    parser.add_argument("--out_parquet", help="Write the MultiIndex dataframe to this output file so it can be used by other scripts (parquet)", required=False)
    parser.add_argument("--include_raw", action='store_true', help="Include a sheet with the raw data", required=False, default=False)
    
    parser.add_argument("--version", action="version", version="0.0.1")

    return parser.parse_args()

    
def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    jc = JoinAndCompare()
    
    dataframes = []
    
    # Some dataframes are optional
    if args.tfx_nomenclature:
        dataframes.append(jc.get_nomenclature_df(NomenclatureTools.TFX.value, args.tfx_nomenclature))
    if args.cgd_nomenclature:
        dataframes.append(jc.get_nomenclature_df(NomenclatureTools.CGD.value, args.cgd_nomenclature))
    if args.variant_validator_nomenclature:
        dataframes.append(jc.get_nomenclature_df(NomenclatureTools.VARIANT_VALIDATOR.value, args.variant_validator_nomenclature))
    if args.mutalyzer_nomenclature:
        dataframes.append(jc.get_nomenclature_df(NomenclatureTools.MUTALYZER.value, args.cgd_nomenclature))
    if args.snpeff_nomenclature:
        dataframes.append(jc.get_nomenclature_df(NomenclatureTools.SNPEFF.value, args.snpeff_nomenclature))
            
    # In the future i will make this optional but for now they're always being generated 
    dataframes.append(jc.get_nomenclature_df(NomenclatureTools.ANNOVAR.value, args.annovar_nomenclature, '.annovar'))    
    
    dataframes.append(jc.get_nomenclature_df(NomenclatureTools.VEP_REFSEQ.value, args.vep_refseq_nomenclautre, 'vep.refseq.'))    
    dataframes.append(jc.get_nomenclature_df(NomenclatureTools.VEP_HG19.value, args.vep_hg19_nomenclature, 'vep.hg19.'))
    
    # Transcripts with alignment gaps 
    gff_and_uta_exon_gap_info_df = jc.get_info_df(args.gff_and_uta_exon_gap_info, NomenclatureTools.REFERENCE_GAP.value)
    
    # Transcripts that have been reported in CGD
    preferred_transcript_df = jc.get_info_df(args.preferred_transcripts, NomenclatureTools.PREFERRED_TRANSCRIPT.value)
    
    # Compare exon, c., and p. from all datasources. 
    comparison_df = jc.get_comparison_df(dataframes, gff_and_uta_exon_gap_info_df, preferred_transcript_df)
        
    jc._generate_concordance_tables(comparison_df)
    
    # Left join all datasets to the variant list. Write out all rows and all fields
    jc.write(args.out, comparison_df, args.include_raw)
    

if __name__ == '__main__':
    main()
