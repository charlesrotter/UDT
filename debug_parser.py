"""
Debug SPARC Parsing - Find the exact issue
"""

def debug_sparc_line():
    """Debug the exact parsing issue step by step"""
    
    test_line = "     D631-7 10   7.72  0.18  2 59.0  3.0   0.196   0.009  1.22    20.93  0.70   115.04   0.290  0.00  57.7   2.7   1      Tr09,dB01"
    
    print("🔍 DEBUGGING SPARC LINE PARSING")
    print("=" * 60)
    print(f"Line length: {len(test_line)}")
    print(f"Line: '{test_line}'")
    print()
    
    # Show character positions
    print("Character positions (every 10):")
    for i in range(0, len(test_line), 10):
        print(f"{i:3d}: '{test_line[i:i+10]}'")
    print()
    
    # Test each field extraction
    fields_to_test = [
        ('name', 0, 11, "Galaxy Name"),
        ('type', 11, 13, "Hubble Type"),
        ('distance', 13, 19, "Distance"),
        ('distance_err', 19, 24, "Distance Error"),
        ('inc', 26, 30, "Inclination"),
        ('inc_err', 30, 34, "Inc Error"),
        ('luminosity', 34, 41, "L[3.6]"),
        ('luminosity_err', 41, 48, "L[3.6] Error"),
        ('reff', 48, 53, "Effective Radius"),
        ('vflat', 86, 91, "V_flat"),
        ('vflat_err', 91, 96, "V_flat Error"),
        ('quality', 96, 99, "Quality Flag")
    ]
    
    print("Field extraction test:")
    print("-" * 50)
    
    for field_name, start, end, description in fields_to_test:
        try:
            raw_value = test_line[start:end]
            stripped_value = raw_value.strip()
            
            print(f"{field_name:12s} [{start:2d}:{end:2d}]: '{raw_value}' -> '{stripped_value}'")
            
            # Try to convert to appropriate type
            if field_name in ['name']:
                converted = stripped_value
            elif field_name in ['type', 'quality']:
                converted = int(stripped_value) if stripped_value else 0
            else:
                converted = float(stripped_value) if stripped_value else 0.0
                
            print(f"                    Converted: {converted}")
            
        except Exception as e:
            print(f"{field_name:12s} [{start:2d}:{end:2d}]: ERROR - {e}")
        
        print()
    
    # Try the actual parsing function
    print("="*60)
    print("TESTING ACTUAL PARSING FUNCTION:")
    
    try:
        galaxy = parse_sparc_line_debug(test_line)
        if galaxy:
            print("✅ PARSING SUCCESSFUL!")
            for key, value in galaxy.items():
                print(f"  {key}: {value}")
        else:
            print("❌ PARSING RETURNED None")
    except Exception as e:
        print(f"❌ PARSING EXCEPTION: {e}")
        import traceback
        traceback.print_exc()

def parse_sparc_line_debug(line: str) -> dict:
    """Debug version of parsing function with detailed error reporting"""
    
    print(f"Input line length: {len(line)}")
    
    if len(line) < 113:
        print(f"❌ Line too short: {len(line)} < 113")
        return None
        
    try:
        # Use exact byte positions (convert to 0-based indexing)
        galaxy = {}
        
        # Parse each field with error checking
        try:
            galaxy['name'] = line[0:11].strip()
            print(f"✅ name: '{galaxy['name']}'")
        except Exception as e:
            print(f"❌ name error: {e}")
            return None
            
        try:
            type_str = line[11:13].strip()
            galaxy['type'] = int(type_str) if type_str else 0
            print(f"✅ type: {galaxy['type']}")
        except Exception as e:
            print(f"❌ type error: {e}")
            return None
            
        try:
            dist_str = line[13:19].strip()
            galaxy['distance'] = float(dist_str) if dist_str else 0.0
            print(f"✅ distance: {galaxy['distance']}")
        except Exception as e:
            print(f"❌ distance error: {e}")
            return None
            
        try:
            lum_str = line[34:41].strip()
            galaxy['luminosity'] = float(lum_str) if lum_str else 0.0
            print(f"✅ luminosity: {galaxy['luminosity']}")
        except Exception as e:
            print(f"❌ luminosity error: {e}")
            return None
            
        try:
            vflat_str = line[86:91].strip()
            galaxy['vflat'] = float(vflat_str) if vflat_str else 0.0
            print(f"✅ vflat: {galaxy['vflat']}")
        except Exception as e:
            print(f"❌ vflat error: {e}")
            return None
            
        try:
            vflat_err_str = line[91:96].strip()
            galaxy['vflat_err'] = float(vflat_err_str) if vflat_err_str else 0.0
            print(f"✅ vflat_err: {galaxy['vflat_err']}")
        except Exception as e:
            print(f"❌ vflat_err error: {e}")
            return None
        
        # Add remaining required fields
        galaxy['distance_err'] = 0.0
        galaxy['inc'] = 0.0
        galaxy['inc_err'] = 0.0
        galaxy['luminosity_err'] = 0.0
        galaxy['reff'] = 0.0
        galaxy['sbeff'] = 0.0
        galaxy['rdisk'] = 0.0
        galaxy['sbdisk'] = 0.0
        galaxy['mhi'] = 0.0
        galaxy['rhi'] = 0.0
        galaxy['quality'] = 2
        
        # Clean up invalid values
        if galaxy['vflat'] <= 0:
            galaxy['vflat'] = float('nan')
        if galaxy['vflat_err'] <= 0 or galaxy['vflat_err'] > 900:
            galaxy['vflat_err'] = float('nan')
            
        return galaxy
        
    except Exception as e:
        print(f"❌ Overall parsing error: {e}")
        import traceback
        traceback.print_exc()
        return None

# Alternative: Try whitespace-based parsing
def parse_sparc_line_whitespace(line: str) -> dict:
    """Alternative parsing using whitespace separation"""
    
    print("\n🔄 TRYING WHITESPACE-BASED PARSING:")
    
    parts = line.split()
    print(f"Split into {len(parts)} parts:")
    for i, part in enumerate(parts):
        print(f"  {i:2d}: '{part}'")
    
    if len(parts) < 15:
        print("❌ Not enough parts")
        return None
    
    try:
        galaxy = {
            'name': parts[0],
            'type': int(parts[1]),
            'distance': float(parts[2]),
            'distance_err': float(parts[3]),
            'inc': float(parts[5]),
            'inc_err': float(parts[6]),
            'luminosity': float(parts[7]),
            'luminosity_err': float(parts[8]),
            'reff': float(parts[9]),
            'sbeff': float(parts[10]),
            'rdisk': float(parts[11]),
            'sbdisk': float(parts[12]),
            'mhi': float(parts[13]),
            'rhi': float(parts[14]),
            'vflat': float(parts[15]),
            'vflat_err': float(parts[16]),
            'quality': int(parts[17])
        }
        
        print("✅ WHITESPACE PARSING SUCCESSFUL!")
        print(f"  Galaxy: {galaxy['name']}")
        print(f"  Luminosity: {galaxy['luminosity']}")
        print(f"  Vflat: {galaxy['vflat']} ± {galaxy['vflat_err']}")
        
        return galaxy
        
    except Exception as e:
        print(f"❌ Whitespace parsing error: {e}")
        return None

if __name__ == "__main__":
    debug_sparc_line()
    
    # Also try whitespace parsing
    test_line = "     D631-7 10   7.72  0.18  2 59.0  3.0   0.196   0.009  1.22    20.93  0.70   115.04   0.290  0.00  57.7   2.7   1      Tr09,dB01"
    parse_sparc_line_whitespace(test_line)