from sparc_data_parser import SPARCDatabaseParser

# Test the parser
parser = SPARCDatabaseParser()

# Test single line
test_line = "     D631-7 10   7.72  0.18  2 59.0  3.0   0.196   0.009  1.22    20.93  0.70   115.04   0.290  0.00  57.7   2.7   1      Tr09,dB01"
galaxy = parser._parse_sparc_line(test_line)

if galaxy:
    print("✅ PARSING WORKS!")
    print(f"Galaxy: {galaxy['name']}")
    print(f"Vflat: {galaxy['vflat']} km/s")
else:
    print("❌ PARSING FAILED")