def debug_sparc_format():
    with open("data/sparc_database/SPARC_Lelli2016c.mrt", 'r') as f:
        lines = f.readlines()
    
    print("🔍 DEBUGGING SPARC FILE FORMAT")
    
    # Find first data line
    for i, line in enumerate(lines):
        if 'CamB' in line or 'DDO' in line:
            print(f"Line {i}: {line}")
            print(f"Length: {len(line)}")
            parts = line.split()
            print(f"Parts: {parts}")
            break

# Run this first
debug_sparc_format()