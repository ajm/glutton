import os

class LocalInfo(object) :
    def __init__(self, options) :
        self.dbdir = options['dbdir']
        self.verbose = options['verbose']

        self.caches = self._get_caches()

    def _get_caches(self) :
        tmp = {}

        if not os.path.exists(self.dbdir) :
            return tmp

        for releasedir in os.listdir(self.dbdir) :
            try :
                release = int(releasedir)
            except ValueError :
                continue

            for speciesdir in os.listdir(os.path.join(self.dbdir, releasedir)) :
                done = os.path.exists(os.path.join(self.dbdir, releasedir, speciesdir, 'done'))

                if speciesdir not in tmp :
                    tmp[speciesdir] = []

                tmp[speciesdir].append((release, done))

        return tmp

    def print_species_table(self) :
        
        if len(self.caches) == 0 :
            print "no databases found in %s" % self.dbdir
            return

        n_len = len(sorted(self.caches, key=len, reverse=True)[0]) + 2
        r_len = len("Release") + 2
        c_len = len("Complete") + 2
    
        n_len = 31
        r_len = 20
        c_len = 10

        print ""
        print "Name".rjust(n_len) + "Release".rjust(r_len) + "Complete".rjust(c_len)
        print "-" * (n_len + r_len + c_len)

        for species in sorted(self.caches) :
            for release, done in sorted(self.caches[species]) :
                print species.rjust(n_len) + str(release).rjust(r_len) + ("Y" if done else "N").rjust(c_len)

    def is_valid_species(self, species) :
        return species in self.caches

    def is_valid_release(self, species, release) :
        if species not in self.caches :
            return False

        for rel, done in self.caches[species] :
            if rel == release :
                return True

        return False

    def get_latest_release(self, species) :
        try :
            return sorted(self.caches[species])[-1][0]
        except :
            return -1

