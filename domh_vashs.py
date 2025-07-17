CREATE TABLE acmg_criteria (
    criterion_id SERIAL PRIMARY KEY,
    name VARCHAR(10) NOT NULL,  -- 'PS1', 'PM5', κλπ.
    description TEXT,
    type VARCHAR(10)  -- 'Pathogenic', 'Benign'
);

INSERT INTO acmg_criteria (name, description, type) VALUES
('PS1', 'Same amino acid change as known pathogenic variant', 'Pathogenic'),
('PM5', 'Different amino acid change at known pathogenic position', 'Pathogenic'),
('PP5', 'Pathogenic classification from trusted source', 'Pathogenic'),
('BP6', 'Benign classification from trusted source', 'Benign');