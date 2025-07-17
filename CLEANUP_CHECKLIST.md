# Pre-Publication Cleanup Checklist

Use this checklist before publishing or sharing the UDT repository.

## 🧹 Cleanup Steps

### 1. Experimental Code
- [ ] Review all files in `exploration/` directory
- [ ] Move successful experiments to main codebase with proper names
- [ ] Document important dead ends in main documentation
- [ ] Archive exploration directory: `tar -czf exploration_archive_YYYY-MM-DD.tar.gz exploration/`
- [ ] Remove exploration directory: `rm -rf exploration/`

### 2. File Cleanup
- [ ] Remove all `experimental_*.py` files
- [ ] Remove all `draft_*.py` files  
- [ ] Remove all `test_idea_*.py` files
- [ ] Remove any `*.tmp` or `*.bak` files
- [ ] Clean up any large data files not needed

### 3. Git History Cleanup (Optional - Destructive!)
```bash
# Create backup first!
git checkout -b pre-cleanup-backup

# Remove specific files from history
git filter-branch --tree-filter 'rm -rf exploration/' HEAD

# Or use git filter-repo (recommended)
pip install git-filter-repo
git filter-repo --path exploration/ --invert-paths
```

### 4. Results Cleanup
- [ ] Remove failed analysis results
- [ ] Consolidate duplicate results
- [ ] Ensure all results have proper documentation
- [ ] Remove any results with errors

### 5. Documentation Review
- [ ] Update README.md to remove references to experimental work
- [ ] Ensure CLAUDE.md is current
- [ ] Remove any TODO comments that are obsolete
- [ ] Check all links and references

### 6. Code Quality
- [ ] Run linting on all Python files
- [ ] Ensure consistent naming conventions
- [ ] Add missing docstrings
- [ ] Remove commented-out code blocks

### 7. Data Files
- [ ] Verify all data files are necessary
- [ ] Add large data files to .gitignore
- [ ] Consider using Git LFS for large files
- [ ] Document data sources clearly

### 8. Final Checks
- [ ] Run all main analysis scripts to ensure they work
- [ ] Check that installation instructions are complete
- [ ] Verify all dependencies in requirements.txt
- [ ] Test on clean environment

## 🚀 Publication Branch

Create clean publication branch:
```bash
git checkout -b publication
# Do cleanup
git add .
git commit -m "Prepare for publication"
```

## 📝 What to Keep

Always keep:
- Core implementation files
- Validated analysis scripts  
- Documentation of the journey (what worked AND what didn't)
- Test files that validate the theory
- Clear data providence

## ⚠️ Before Destructive Operations

1. **Always create backup branch**: `git checkout -b backup-before-cleanup`
2. **Push to private remote**: `git push origin backup-before-cleanup`
3. **Local archive**: `tar -czf udt-complete-backup.tar.gz .`

---

*Last cleanup performed: Never (repository still in development)*