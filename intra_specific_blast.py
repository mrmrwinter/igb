def file_translate(input_fasta):
    from Bio import SeqIO
    records = list(SeqIO.parse(input_fasta, 'fasta'))
    translated = []
    for record in records:
        translated.append(record.seq.translate())
    SeqIO.write(translated, input_fasta+'.aa','fasta')
    return input_fasta+'.aa'

class IntraSpecificBlastCommandLine:
    
    def __init__(self, cdss, makeblastdb='makeblastdb', cd_hit_est='cd-hit-est', **kwargs):
        
        cd_hit_cline = cd_hit_est+' -i '+cdss+' -o '+cdss+'.98 -c 0.98 -g 1'
        
        self.cd_hit_est_cline = cd_hit_cline
        
        db_cline = makeblastdb+' -in '+cdss+'.98 -dbtype nucl'
        self.makeblastdb_cline = db_cline

        from Bio.Blast.Applications import NcbiblastnCommandline
        
        kwargs['query'] = cdss+'.98'
        kwargs['db'] = cdss+'.98'
        kwargs['outfmt'] = 5
        
        cline = NcbiblastnCommandline(**kwargs)
        
        self.blastn_cline = cline
        
        
        self.results_handle = None
        
    def execute(self):
        
        import subprocess as sub
        
        sub.call(self.cd_hit_est_cline, shell=True)
        
        import time
        time.sleep(30)
        
        sub.call(self.makeblastdb_cline, shell=True)
        time.sleep(30)
        
        self.results = self.blastn_cline()[0]
        
    def calculate_percents(self):
        
        from Bio.Blast import NCBIXML
        import io

        string_results = io.BytesIO(self.results)
        results = NCBIXML.parse(string_results)
        
        percents = []
        
        self.percent_per_hit = {}
        
        for query in results:
            query_id = query.query.split()[0]
            query_length = query.query_length
            percent = None
            for aln in query.alignments:
                match_id = aln.title.split('|')[-1].split()[1]
                for hsp in aln.hsps:
                    target_length = len(str(hsp.sbjct).replace('-',''))
                    prop_target_of_query = float(target_length)/int(query_length)
                    identities = hsp.identities
                    if not percent and not match_id == query_id and prop_target_of_query > 0.7:
                        percent = 100*(float(identities)/int(target_length))
                        percents.append(percent)
                        self.percent_per_hit[query_id] = {'hit': match_id,
                                                          'percent_ident': percent,
                                                          'prop_query': prop_target_of_query}
        self.percent_per_hit            
        return percents
    
    def write_blast_results(self, lower_cutoff = 0, upper_cutoff = 100):
        report = ''
        report += ('matches on other contigs than query with %s <= identity <= %s\n'%(str(lower_cutoff), str(upper_cutoff))).upper()
        report += '------------------------------------------------------------\n'
        report += '\n'
        A = self.percent_per_hit
        maxlength = max([len(q) for q in A.keys()])
        report += 'query'.ljust(maxlength+1)+'target'.ljust(maxlength+1)+'%identity '+'prop_lenght(match/query)\n'
        report += '\n'
        count = 0
        for query in A:
            q=query
            if (A[q]['percent_ident'] >= lower_cutoff and
                A[q]['percent_ident'] <= upper_cutoff):
                report += (q.ljust(maxlength+1)+A[q]['hit'].ljust(maxlength+1)+
                           str(A[q]['percent_ident'])[:6].ljust(10)+
                           str(A[q]['prop_query'])[:6].ljust(10)+'\n')
                count += 1
        report += '\n'
        report += str(count)+'\n'
        report += 'Total matches: ' + str(len(A.keys()))
        return report

