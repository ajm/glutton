import datetime
from sys import stderr
from xml.sax import ContentHandler, make_parser
from xml.sax.saxutils import XMLGenerator

from glutton.genefamily import ensembl_to_glutton_internal, Gene, GeneFamily


timefmt = '%Y-%m-%dT%H:%M:%S'


class GluttonTypeNotSupportedError(Exception) :
    pass

class GluttonXMLManifestHandler(ContentHandler) :
    def __init__(self) :
        self.metadata = {}

    def startElement(self, name, attrs) :
        if name == 'manifest' :
            return

        self.current_field = name
        self.current_type = attrs['type']

    def characters(self, content) :
        self.metadata[self.current_field] = self._convert(content, self.current_type)

    def _convert(self, content, t) :
        global timefmt

        if t == 'str' or t == 'unicode':
            return content
        elif t == 'int' :
            return int(content)
        elif t == 'float' :
            return float(content)
        elif t == 'bool' :
            return content == 'True'
        elif t == 'datetime' :
            return datetime.datetime.strptime(content, timefmt)

        raise GluttonTypeNotSupportedError("unknown type: %s" % t)

class GluttonXMLManifest(object) :
    
    def read(self, f) :
        self.parser = make_parser()
        self.parser.setContentHandler(GluttonXMLManifestHandler())
        
        self.parser.parse(f)

        return self.parser.getContentHandler().metadata

    def write(self, f, metadata) :
        global timefmt
        
        writer = XMLGenerator(f, 'utf-8')

        writer.startDocument()
        writer.startElement('manifest', {})

        for m in metadata :
            t = type(metadata[m]).__name__
            writer.startElement(m, { 'type' : t })

            # everything else is a builtin type, but datetime is 
            # good to have as well
            if t == 'datetime' :
                content = metadata[m].strftime(timefmt)
            else :
                content = str(metadata[m])

            writer.characters(content)
            writer.endElement(m)

        writer.endElement('manifest')
        writer.endDocument()

        f.flush()

class GluttonXMLMappingHandler(ContentHandler) :
    def __init__(self) :
        self.mapping = {}

    def startElement(self, name, attrs) :
        if name == 'mapping' :
            self.current_key = attrs['id']

    def characters(self, content) :
        self.current_value = content

    def endElement(self, name) :
        if name == 'mapping' :
            self.mapping[self.current_key] = self.current_value

class GluttonXMLMapping(object) :

    def read(self, f) :
        self.parser = make_parser()
        self.parser.setContentHandler(GluttonXMLMappingHandler())

        self.parser.parse(f)

        return self.parser.getContentHandler().mapping

    def write(self, f, mapping) :
        writer = XMLGenerator(f, 'utf-8')

        writer.startDocument()
        writer.startElement('mappings', {})

        for m in mapping :
            writer.startElement('mapping', { 'id' : m })
            writer.characters(mapping[m])
            writer.endElement('mapping')

        writer.endElement('mappings')
        writer.endDocument()

        f.flush()

class GluttonXMLDataHandler(ContentHandler):
    def __init__(self) :
        self.families = {}

        self.current_family = None
        self.current_gene = None

    def get_data(self) :
        return self.families

    def startElement(self, name, attrs):
        if name == 'genefamiles' :
            pass

        elif name == 'genefamily' :
            self.current_family = GeneFamily(id=attrs['id'])

        elif name == 'gene' :
            self.current_gene = Gene(attrs['name'], id=attrs['id'])

    def endElement(self, name) :
        if name == 'genefamiles' :
            pass

        elif name == 'genefamily' :
            self.families[self.current_family.id] = self.current_family

        elif name == 'gene' :
            self.current_family.append(self.current_gene)

    def characters(self, content) :
        self.current_gene.sequence = content

class GluttonXMLData(object) :
    
    def read(self, f) :
        self.parser = make_parser()
        self.parser.setContentHandler(GluttonXMLDataHandler())
        
        self.parser.parse(f)

        ch = self.parser.getContentHandler()
        
        return ch.get_data()

    def write(self, f, families) :
        writer = XMLGenerator(f, 'utf-8')
        
        writer.startDocument()
        writer.startElement('genefamilies', {})

        for famid in families :
            writer.startElement('genefamily', { 'id' : famid })

            for p in families[famid] :
                writer.startElement('gene', { 'id' : p.id, 'name' : p.name })
                writer.characters(p.sequence)
                writer.endElement('gene')

            writer.endElement('genefamily')

        writer.endElement('genefamilies')
        writer.endDocument()

        f.flush()


if __name__ == '__main__' :
    # XML manifest
    import sys
    import glutton
    import datetime

    metadata = {
            'glutton-version'   : glutton.__version__,
            'program-name'      : 'PRANK',
            'program-version'   : 'v1.1',
            'species-name'      : 'tribolium_castaneum',
            'species-release'   : 23,
            'download-time'     : datetime.datetime.today(),
            'data-file'         : 'tribolium_castaneum_23.xml',
            'mapping-file'      : 'tc23_map.xml',
            'is-complete'       : True
        }

    gp = GluttonXMLManifest()
    gp.write(open('manifest.xml', 'w'), metadata)

    print gp.read(open('manifest.xml', 'r'))


if __name__ == '__main__x' :
    # XML data
    from glutton.ensembl import EnsemblDownloader

    e = EnsemblDownloader()
    fam = e.download('tribolium_castaneum', 23)
    print "downloaded %d gene families" % len(fam)

    gp = GluttonXMLData()
    fname = 'tribolium_castaneum_test.xml'

    gp.write(open(fname, 'w'), ensembl_to_glutton_internal(fam))
    print "written %s" % fname
    
    fam2 = gp.read(open(fname, 'r'))
    print "read %d families from %s" % (len(fam2), fname)

