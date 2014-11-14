from sys import exit

from glutton.utils import get_log

from glutton.ensembl_biomart import download_database_biomart, get_all_species_biomart, get_latest_release_biomart
from glutton.ensembl_pycogent import download_database_pycogent, get_all_species_pycogent
from glutton.ensembl_sql import download_database_sql, get_all_species_sql, get_latest_release_sql, ReleaseNotFoundError, SQLQueryError, NoResultsError


ENSEMBL_METHODS = ('biomart', 'sql', 'pycogent')
DEFAULT_METHOD = 'sql'

def set_ensembl_download_method(method) :
    global DEFAULT_METHOD, ENSEMBL_METHODS

    assert method in ENSEMBL_METHODS

    DEFAULT_METHOD = method

def get_ensembl_download_method() :
    global DEFAULT_METHOD

    return DEFAULT_METHOD

class EnsemblDownloadError(Exception) :
    pass

class EnsemblDownloader(object) :
    def __init__(self) :
        self.log = get_log()
        self.method = get_ensembl_download_method()

        self.log.info("ensembl download method is set to: %s" % self.method)

    # XXX   i don't want the calling code to do this anymore, I want it to be handled by
    #       the download code if Release is None

    def get_latest_release(self, species) :
        if self.method in ('sql', 'pycogent') :
            return get_latest_release_sql(species)

        elif self.method == 'biomart' :
            return get_latest_release_biomart(species)

        else :
            self.log.fatal("%s not a valid download method" % method)
            exit(1)

    def get_all_species(self, db="", suppress=None) :

        if self.method == 'biomart' :
            tmp = get_all_species_biomart(db, suppress)
        
        elif self.method == 'sql' :
            tmp = get_all_species_sql(db, suppress)
        
        elif self.method == 'pycogent' :
            tmp = get_all_species_pycogent(db, suppress)
        
        else :
            self.log.fatal("%s not a valid download method" % method)
            exit(1)


        return sorted(tmp)

    # returns peptide sequences + homologies
    def download(self, species, release=None, nucleotide=False, slow=False) :
        try :
            if self.method == 'biomart' :
                return download_database_biomart(species, release, nucleotide)

            elif self.method == 'sql' :
                return download_database_sql(species, release, nucleotide)

            elif self.method == 'pycogent' :
                return download_database_pycogent(species, release, nucleotide)

            else :
                raise EnsemblDownloadError("%s is not a valid download method" % self.method)

        except ReleaseNotFoundError, rnfe :
            raise EnsemblDownloadError(rnfe.message)

        except SQLQueryError, sql :
            raise EnsemblDownloadError(' '.join(sql.message.replace('\\n', '').split()))

        except NoResultsError :
            raise EnsemblDownloadError("no results returned from database")

#        except UnsupportedError, ue :
#            raise EnsemblDownloadError(ue.message)

        except Exception, e :
            raise EnsemblDownloadError(e.message + " (" + self.method + ")")

