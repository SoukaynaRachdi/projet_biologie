const analyzeBtn = document.getElementById('analyze-btn');
const dnaInput = document.getElementById('dna-input');
const resultBox = document.getElementById('result-box');
const fastaFile = document.getElementById('fasta-file');

fastaFile.addEventListener('change', (e) => {
    const file = e.target.files[0];
    const reader = new FileReader();
    reader.onload = function() {
        const lines = reader.result.split('\n');
        let seq = '';
        for (let line of lines) {
            if (!line.startsWith('>')) seq += line.trim();
        }
        dnaInput.value = seq;
    };
    if(file) reader.readAsText(file);
});

analyzeBtn.addEventListener('click', async () => {
    const dna = dnaInput.value.trim().toUpperCase();
    if (!dna) {
        alert("Please enter a DNA sequence.");
        return;
    }

    // Appel à Flask backend
    const response = await fetch('/analyze', {
        method: 'POST',
        headers: {'Content-Type':'application/json'},
        body: JSON.stringify({dna})
    });

    const data = await response.json();
    resultBox.value = data.result;
});