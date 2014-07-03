  use XML::SAX::ParserFactory;
  use XML::Validator::Schema;

  #
  # create a new validator object, using foo.xsd
  #
  $validator = XML::Validator::Schema->new(file => 'schemas/nb.xsd');

  #
  # create a SAX parser and assign the validator as a Handler
  #
  $parser = XML::SAX::ParserFactory->parser(Handler => $validator);

  #
  # validate foo.xml against foo.xsd
  #
  eval { $parser->parse_uri('SBM/sbm.nb') };
  die "File failed validation: $@" if $@;
