#include <iostream>

#include <QCoreApplication>
#include <QCommandLineOption>
#include <QCommandLineParser>
#include <QDebug>
#include <QDirIterator>
#include <QFile>
#include <QTextCodec>

using QCLO = QCommandLineOption;

QString process(QString line, QString var, int lineNO, int& nReplaces, int& nWarnings, int& nErrors) {
    if(line.count(var, Qt::CaseInsensitive) >= 2) {
        nErrors++;
        qDebug() << "\nERROR: Variable" << var << "kommt mehr als einmal in dieser Zeile vor. Append >>>>>> to the line";
        qDebug() << lineNO << ":" << line << "\n";
        return ">>>>>>" + line;
    }

    int start = line.indexOf(var, 0, Qt::CaseInsensitive);
    if(start == -1) return line;

    start += var.size();

    // Das erste non-whitespace Zeichen nach der Variable muss einen oeffnende runde Klammer sein.
    // Ansonsten ist es irgendein anderer Ausdruck den wir nicht anedern wollen.
    bool ja = false;
    for(int i=start; i < line.size(); ++i) {
        if(line[i] == ' ') continue;
        if(line[i] == '(') {
            ja = true;
        }
        break;
    }
    if(not ja) return line;

    int LHSstart = 0;
    int LHSend = 0;
    int RHSstart = 0;
    int RHSend = 0;

    int tiefe = 0;
    for(int i=start; i < line.size(); ++i) {
        if(line[i] == '(') {
            if(tiefe == 0) {
                LHSstart = i+1;
            }

            tiefe++;

        } else if(line[i] == ')') {
            tiefe--;

            if(tiefe == 0) {
                RHSend = i;
                break;
            }
        } else if(tiefe == 1 and line[i] == ',') {
            LHSend = i;
            RHSstart = i+1;
        }
    }

    if(LHSstart==0 or LHSend==0 or RHSstart==0 or RHSend==0) {
        nErrors++;
        qDebug() << "\nError: Klammern Pattern nicht gefunden. Append >>>>>> to the line";
        qDebug() << lineNO << ":" << line << "\n";
        return ">>>>>>" + line;
    }

    nReplaces++;

    QString newline = line.left(LHSstart) + line.mid(RHSstart, RHSend-RHSstart) + "," + line.mid(LHSstart, LHSend-LHSstart) + line.right(line.size()-RHSend);
//    qDebug() << "old: " << line;
//    qDebug() << "new: " << newline;
    return newline;
}

int main(int argc, char *argv[]) {
    try {
        QCoreApplication a(argc, argv);
        qDebug() << QTextCodec::codecForLocale()->name();

        QCLO inputDirOption(  "input-dir",    "Input directory",                                                     "path");
        QCLO outputDirOption( "output-dir",   "Output directory",                                                    "path");
        QCLO variableOption(  "variable",     "Variable wich idx to replace (Option can be given multiple times)",   "name");

        QCommandLineParser parser;
        parser.addOption(inputDirOption);
        parser.addOption(outputDirOption);
        parser.addOption(variableOption);

        parser.process(a);

        if(not parser.isSet(inputDirOption)) throw "Option --input-dir must be given";
        if(not parser.isSet(outputDirOption)) throw "Option --output-dir must be given";
        if(not parser.isSet(variableOption)) throw "Option --variable must be given";

        QString inputDir = parser.value(inputDirOption);
        QString outputDir = parser.value(outputDirOption);
        QStringList variables = parser.values(variableOption);

        qDebug() << "Input directory: " << inputDir;
        qDebug() << "Output directory: " << outputDir;
        qDebug() << "Variables: " << variables;



        int totalFiles = 0;
        int totalLines = 0;
        int totalWarning = 0;
        int totalErrors = 0;
        int totalReplaces = 0;

        QDirIterator dirIter(inputDir);
        while(dirIter.hasNext()) {
            QString nextFile = dirIter.next();
            QString fileName = QFileInfo(nextFile).fileName();
            qDebug() << "process:" << fileName;

            QFile file(nextFile);
            if(not file.open(QFile::ReadOnly | QFile::Text)) {
                totalErrors++;
                qDebug() << "ERROR: Can not open file for read" << nextFile;
                continue;
            }

            QFile outFile(outputDir+"/"+fileName);
            if(not outFile.open(QFile::WriteOnly | QFile::Text)) {
                totalErrors++;
                qDebug() << "ERROR: Can not open file for write" << outputDir+"/"+fileName;
                continue;
            }
            QTextStream out(&outFile);
            // $ file mod_xnl4v5.f90
            // mod_xnl4v5.f90: ISO-8859 text
            out.setCodec(QTextCodec::codecForName("ISO 8859-1"));


            int lineNo = 0;
            int nReplaces = 0;
            QTextStream in(&file);
            in.setCodec(QTextCodec::codecForName("ISO 8859-1"));
            while (not in.atEnd()) {
                QString line = in.readLine();
                lineNo++;
                for(QString var : variables) {
                   line = process(line, var, lineNo, nReplaces, totalWarning, totalErrors);
                }
                out << line << "\n";
            }

            if(nReplaces > 0) qDebug() << nReplaces << "Ersetzungen";

            totalFiles++;
            totalLines += lineNo;
            totalReplaces += nReplaces;
        }

        qDebug() << "\n";
        qDebug() << "Processed" << totalFiles << "file with" << totalLines << "lines";
        qDebug() << "Replaces:" << totalReplaces;
        qDebug() << "Errors  :" << totalErrors;
        qDebug() << "Warnings:" << totalWarning;

        return EXIT_SUCCESS;
    } catch(const char* ex) {
        std::cerr << ex << "\n";
    }
    return EXIT_FAILURE;
}
